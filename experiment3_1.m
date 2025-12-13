%---这是近焦点转地心惯性坐标系---%



% 轨道参数
a = 6700;
e = 0.5;
i = deg2rad(30);
Omega = deg2rad(30);
omega = deg2rad(30);

% 时间向量
t = 0:100000:100000000000;

% 平均运动
n = sqrt(1 / (a^3));

% 计算偏近点角
E = zeros(size(t));
E(1) = 0;
for k = 2:length(t)
    M = n * t(k);
    E(k) = M + e * sin(E(k-1));
end

% 真近点角
theta = 2 * atan(sqrt((1 + e) / (1 - e)) * tan(E / 2));

% 轨道半径
r = a * (1 - e^2) ./ (1 + e * cos(theta));

% 近地点坐标系位置
r_pqw = [r .* cos(theta); r .* sin(theta); zeros(size(r))];

% 旋转矩阵
R = [cos(Omega)*cos(omega) - sin(Omega)*sin(omega)*cos(i), ...
     -cos(Omega)*sin(omega) - sin(Omega)*cos(omega)*cos(i), ...
      sin(Omega)*sin(i);
     sin(Omega)*cos(omega) + cos(Omega)*sin(omega)*cos(i), ...
     -sin(Omega)*sin(omega) + cos(Omega)*cos(omega)*cos(i), ...
     -cos(Omega)*sin(i);
     sin(omega)*sin(i), ...
     cos(omega)*sin(i), ...
     cos(i)];

% 转换到惯性坐标系
r_inertial = zeros(length(t), 3);
for k = 1:length(t)
    r_inertial(k, :) = R * r_pqw(:, k);
end

% 可视化
figure('Position', [100, 100, 800, 600]);
hold on;
plot3(r_inertial(:,1), r_inertial(:,2), r_inertial(:,3), 'b-', 'LineWidth', 2);
grid on;
xlabel('X (km)');
ylabel('Y (km)');
zlabel('Z (km)');
title('卫星轨道 (地心惯性坐标系)');
axis equal;
box on;

plot3(0, 0, 0, 'ro', 'MarkerSize', 8, 'LineWidth', 2);
text(0, 0, 0, ' 地球中心', 'FontSize', 10, 'VerticalAlignment', 'bottom');

plot3(r_inertial(1,1), r_inertial(1,2), r_inertial(1,3), 'go', 'MarkerSize', 8, 'LineWidth', 2);
text(r_inertial(1,1), r_inertial(1,2), r_inertial(1,3), ' 起点', 'FontSize', 9, 'VerticalAlignment', 'bottom');

view(3);
camlight('left');
lighting gouraud;
hold off;
