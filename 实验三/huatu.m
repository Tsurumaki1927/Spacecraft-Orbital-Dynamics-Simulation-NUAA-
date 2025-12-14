
%% 数据和常数
load guidaogenshu.mat;

% 地球引力常数
mu = 398600.4415; % km^3/s^2

%% 双曲线限制条件
% 双曲线轨道参数
a_abs = abs(a);  % 双曲线半长轴取绝对值
p = a_abs * (e^2 - 1);  % 半通径

% 确定真近点角范围（避免渐近线处的无穷大）
theta_max = acos(-1/e) * 0.95;  % 接近渐近线但不达到
theta = linspace(-theta_max, theta_max, 1000);

% 双曲线轨道方程：r = p / (1 + e*cosθ)
r_orbit = p ./ (1 + e * cos(theta));

%% 数值求解点
% 双曲线平均运动
n = sqrt(mu / a_abs^3);

% 计算双曲偏近点角F
F = zeros(size(theta));
for k = 1:length(theta)
    % 从真近点角计算双曲偏近点角
    if e > 1
        F(k) = 2 * atanh(sqrt((e-1)/(e+1)) * tan(theta(k)/2));
    end
end

% 计算平近点角M（双曲线开普勒方程）
M = e * sinh(F) - F;

% 计算时间（如果需要）
t = M / n;


% 近地点坐标系中的位置
x_pqw = r_orbit .* cos(theta);
y_pqw = r_orbit .* sin(theta);
z_pqw = zeros(size(theta));

% 旋转矩阵
R = [cos(RAAN)*cos(w) - sin(RAAN)*sin(w)*cos(i), ...
    -cos(RAAN)*sin(w) - sin(RAAN)*cos(w)*cos(i), ...
    sin(RAAN)*sin(i);
    sin(RAAN)*cos(w) + cos(RAAN)*sin(w)*cos(i), ...
    -sin(RAAN)*sin(w) + cos(RAAN)*cos(w)*cos(i), ...
    -cos(RAAN)*sin(i);
    sin(w)*sin(i), ...
    cos(w)*sin(i), ...
    cos(i)];

% 转换到惯性坐标系
r_inertial = zeros(length(theta), 3);
for k = 1:length(theta)
    r_inertial(k, :) = R * [x_pqw(k); y_pqw(k); z_pqw(k)];
end

%% 绘图
figure('Position', [100, 100, 800, 600]);

hold on;
grid on;
box on;

% 绘制轨道
plot3(r_inertial(:,1), r_inertial(:,2), r_inertial(:,3), 'b-', 'LineWidth', 2);
view(3);

% 绘制渐近线


% 标记地心
plot3(0, 0, 0, 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r');

xlabel('X (km)');
ylabel('Y (km)');
zlabel('Z (km)');
title('双曲线轨道 (地心惯性坐标系)');

% 设置坐标轴范围，确保能看到完整的轨道
max_range = max(abs(r_inertial(:))) * 1.2;
xlim([-max_range, max_range]);
ylim([-max_range, max_range]);
zlim([-max_range, max_range]);

axis equal;
view(3);

% 添加图例和说明
legend('双曲线轨道', '地球中心', 'Location', 'best');

