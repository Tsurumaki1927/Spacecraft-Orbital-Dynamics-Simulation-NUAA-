clear; clc; format long g;

%% 1. 物理常数与观测数据
const.mu = 398600.4418;  
const.Re = 6378.137;     
const.f  = 1/298.257223563; 

% 观测站: 北纬60度, 高度0.5km
lat_deg = 60.0;
alt_km  = 0.5;

% 观测数据: 时间(min), LST(deg), RA(deg), Dec(deg)
obs_data = [
    0.0,   150.000,   157.783,  24.2403;
    5.0,   151.253,   159.221,  27.2993;
    10.0,  152.507,   160.526,  29.8982
];

% 数据预处理
t_sec = obs_data(:, 1) * 60;
tau1 = t_sec(1) - t_sec(2); 
tau3 = t_sec(3) - t_sec(2); 
tau  = t_sec(3) - t_sec(1);

LST_rad = deg2rad(obs_data(:, 2)); 
Ra_rad  = deg2rad(obs_data(:, 3)); 
Dec_rad = deg2rad(obs_data(:, 4));

fprintf('tau1=%.0f, tau3=%.0f, tau=%.0f\n', tau1, tau3, tau);

%% 2. 计算测站矢量 R 与 视线矢量 L
R = zeros(3, 3);
L = zeros(3, 3);

sin_phi = sin(deg2rad(lat_deg));
cos_phi = cos(deg2rad(lat_deg));
e2 = 2*const.f - const.f^2;
N = const.Re / sqrt(1 - e2 * sin_phi^2);

R_xy = (N + alt_km) * cos_phi;
R_z  = (N * (1 - e2) + alt_km) * sin_phi;

for i = 1:3
    R(1, i) = R_xy * cos(LST_rad(i));
    R(2, i) = R_xy * sin(LST_rad(i));
    R(3, i) = R_z;
    
    cl = cos(Dec_rad(i));
    L(1, i) = cl * cos(Ra_rad(i));
    L(2, i) = cl * sin(Ra_rad(i));
    L(3, i) = sin(Dec_rad(i));
end

%% 3. 实验内容(1): 高斯法初始解
p1 = cross(L(:,2), L(:,3));
p3 = cross(L(:,1), L(:,2));
D0 = dot(L(:,1), p1);

D = zeros(3,3);
for i = 1:3, D(i,1)=dot(R(:,i),p1); D(i,2)=dot(R(:,i),cross(L(:,1),L(:,3))); D(i,3)=dot(R(:,i),p3); end

A = (1/D0) * (-D(1,2)*(tau3/tau) + D(2,2) + D(3,2)*(tau1/tau));
B = (1/(6*D0)) * (D(1,2)*(tau3^2-tau^2)*(tau3/tau) + D(3,2)*(tau^2-tau1^2)*(tau1/tau));

E = dot(R(:,2), L(:,2)); 
R2_sq = dot(R(:,2), R(:,2));

coeffs = [1, 0, -(A^2 + 2*A*E + R2_sq), 0, 0, -2*const.mu*B*(A + E), 0, 0, -const.mu^2*B^2];
roots_r = roots(coeffs);
valid_r = roots_r(imag(roots_r)==0 & real(roots_r)>const.Re);
[~, idx] = max(valid_r);
r2_mag = valid_r(idx);

% 初始状态计算
u = const.mu / r2_mag^3;
c1 = (tau3/tau) * (1 + u*(tau^2 - tau3^2)/6);
c3 = -(tau1/tau) * (1 + u*(tau^2 - tau1^2)/6); % 修正: c3公式带负号

rho1 = (1/D0)*(-D(1,1) + (1/c1)*D(2,1) - (c3/c1)*D(3,1));
rho2 = A + const.mu * B / r2_mag^3;
rho3 = (1/D0)*(-(c1/c3)*D(1,3) + (1/c3)*D(2,3) - D(3,3));

r1 = R(:,1) + rho1*L(:,1);
r2 = R(:,2) + rho2*L(:,2);
r3 = R(:,3) + rho3*L(:,3);

f1 = 1-0.5*u*tau1^2; g1 = tau1-u*tau1^3/6;
f3 = 1-0.5*u*tau3^2; g3 = tau3-u*tau3^3/6;
v2 = (-f3*r1 + f1*r3) / (f1*g3 - f3*g1);

fprintf('\n(1) 初始解:\n r2 = %.4f km\n v2 = %.4f km/s\n', r2_mag, norm(v2));
fprintf(' r2=(%.4f,%.4f,%.4f) \n',r2);
fprintf(' v2=(%.4f,%.4f,%.4f) \n',v2);

%% 4. 实验内容(2): 迭代改进 (带阻尼)
fprintf('\n(2) 开始迭代...\n');
r2_curr = r2;
lambda = 0.5; % 阻尼因子

for k = 1:30
    r_mag = norm(r2_curr);
    u = const.mu / r_mag^3;
    
    f1 = 1-0.5*u*tau1^2; g1 = tau1-u*tau1^3/6;
    f3 = 1-0.5*u*tau3^2; g3 = tau3-u*tau3^3/6;
    
    det = f1*g3 - f3*g1;
    c1 = g3/det; c3 = -g1/det;
    
    rho1_new = (1/D0)*(-D(1,1) + (1/c1)*D(2,1) - (c3/c1)*D(3,1));
    rho2_new = (1/D0)*(-c1*D(1,2) + D(2,2) - c3*D(3,2)); 
    rho3_new = (1/D0)*(-(c1/c3)*D(1,3) + (1/c3)*D(2,3) - D(3,3));
    
    r1_new = R(:,1) + rho1_new*L(:,1);
    r2_new = R(:,2) + rho2_new*L(:,2);
    r3_new = R(:,3) + rho3_new*L(:,3);
    v2_new = (-f3*r1_new + f1*r3_new) / det;
    
    diff = abs(norm(r2_new) - r_mag);
    fprintf('Iter %2d: r2 = %.4f km, diff = %.4f\n', k, norm(r2_new), diff);
    
    if diff < 1e-8
        r2 = r2_new; v2 = v2_new;
        fprintf('迭代收敛。\n');
        break;
    end
    r2_curr = (1-lambda)*r2_curr + lambda*r2_new; % 阻尼更新
end

fprintf('最终状态:\n r = [%.4f, %.4f, %.4f]\n v = [%.4f, %.4f, %.4f]\n', r2, v2);
%%%这里我选择使用了松弛法求解，避免了发散

%% 5. 实验内容(3): 轨道根数
fprintf('\n(3) 轨道根数:\n');
[a, e, i, RAAN,w , f] = coe_from_sov(r2, v2, const.mu);
save guidaogenshu.mat a e i RAAN w f
fprintf(' a = %.4f km\n e = %.4f\n i = %.4f deg\n', a, e, i);
fprintf(' Ω = %.4f deg\n ω = %.4f deg\n f = %.4f deg\n ', rad2deg(RAAN),rad2deg (w), rad2deg(f));
