close all; clc;clear;tic,

%%%%轨道六根数
mu      = 398600.4418;        % 地球引力常数，单位
R_e     = 6378;               % 地球半径
omega_e = 7.292115e-5;        % 地球自转角速度
h       = 58930;              % 轨道角动量
e       = 0.8;                % 偏心率
incl    = deg2rad(39.67);     % 倾角，直接转弧度
Omega   = deg2rad(130.32);    % 升交点赤经
omega   = deg2rad(42.373);    % 近地点幅角

%%%%提前计算
p = h^2 / mu;
a = p / (1 - e^2);
n = sqrt(mu / a^3);
T = 2*pi / n;

%%%%画图和采样
numPeriods = 3;%三个周期
ptsPerOrbit = 10000;%一万个点提高分辨率，
t = linspace(0, numPeriods*T, numPeriods*ptsPerOrbit + 1).';%计算时间间隔

lat = zeros(size(t));
lon = zeros(size(t));

for k = 1:numel(t)
    M = n * t(k);
    E = keplerE(M, e);
    TA = 2*atan2(sqrt(1+e)*sin(E/2), sqrt(1-e)*cos(E/2));

    [r_eci, ~] = sv_from_coe([h, e, Omega, incl, omega, TA]);

    theta_g = omega_e * t(k);
    R3 = [cos(theta_g)  sin(theta_g) 0;
         -sin(theta_g)  cos(theta_g) 0;
          0             0            1];
    r_ecef = R3 * r_eci;

    r_norm = norm(r_ecef);
    lat(k) = asin(r_ecef(3)/r_norm);
    lon(k) = atan2(r_ecef(2), r_ecef(1));
end

latDeg = rad2deg(lat);
lonDeg = rad2deg(mod(lon + pi, 2*pi) - pi);

figure; hold on; grid on;
colors = lines(numPeriods);
for idx = 1:numPeriods
    seg = (idx-1)*ptsPerOrbit + 1 : idx*ptsPerOrbit + 1;%写一个索引,告诉scatter画到哪个点了
    scatter(lonDeg(seg), latDeg(seg), 12, ...
            'MarkerEdgeColor', colors(idx,:), ...
            'MarkerFaceColor', colors(idx,:), ...
            'MarkerFaceAlpha', 0.6);
end%要绘制散点图，连线图会出错！！！！
xlabel('经度 (°)');
ylabel('纬度 (°)');
title('连续三个周期的星下点轨迹图');
legend('第1圈','第2圈','第3圈','Location','best');
xlim([-180 180]);
ylim([-90 90]);
hold off;
toc


%%%%函数1：开普勒法
function E = keplerE(M, e)
    % %用牛顿法解开普勒方程
    M = mod(M, 2*pi);
    E = M;
    for iter = 1:15%这里迭代了15次
        f  = E - e*sin(E) - M;
        fp = 1 - e*cos(E);
        E  = E - f / fp;
    end
end


%%%%函数2：六根数计算轨道速度和周期
function [r, v] = sv_from_coe(coe) 
    mu = 398600.4418; 
    h = coe(1); 
    e = coe(2); 
    RA = coe(3); 
    incl = coe(4); 
    w = coe(5); 
    TA = coe(6); 
    
    rp = (h^2 / mu) / (1 + e * cos(TA)) * [cos(TA); sin(TA); 0]; 
    vp = (mu / h) * [-sin(TA); e + cos(TA); 0]; 
    
    R3_RA = [cos(RA), sin(RA), 0; -sin(RA), cos(RA), 0; 0, 0, 1]; 
    R1_i = [1, 0, 0; 0, cos(incl), sin(incl); 0, -sin(incl), cos(incl)]; 
    R3_w = [cos(w), sin(w), 0; -sin(w), cos(w), 0; 0, 0, 1]; 
    R = R3_RA * R1_i * R3_w; 
    
    r = R * rp; 
    v = R * vp; 
end
