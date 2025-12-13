%---这是近焦点转地心惯性坐标系在转化回去近焦点---%




a = 6700;
e = 0.5;
i = deg2rad(30);
Omega = deg2rad(30);
omega = deg2rad(30);



t = 0:100000:1000000000;
n = sqrt(1 / (a^3));

M = n * t;
E = M;
for iter = 1:10
    E = M + e * sin(E);
end

theta = 2 * atan(sqrt((1 + e) / (1 - e)) * tan(E / 2));
r = a * (1 - e^2) ./ (1 + e * cos(theta));

r_pqw = [r .* cos(theta); 
         r .* sin(theta); 
         zeros(size(r))];

R = [
    cos(Omega)*cos(omega) - sin(Omega)*sin(omega)*cos(i), ...
    -cos(Omega)*sin(omega) - sin(Omega)*cos(omega)*cos(i), ...
    sin(Omega)*sin(i);
    sin(Omega)*cos(omega) + cos(Omega)*sin(omega)*cos(i), ...
    -sin(Omega)*sin(omega) + cos(Omega)*cos(omega)*cos(i), ...
    -cos(Omega)*sin(i);
    sin(omega)*sin(i), ...
    cos(omega)*sin(i), ...
    cos(i)
];

r_inertial = R * r_pqw;
R_inv = R';
r_pqw_recovered = R_inv * r_inertial;

figure('Position', [100, 100, 800, 600]);
hold on;
plot(r_pqw_recovered(1, :), r_pqw_recovered(2, :), 'b-', 'LineWidth', 2);
grid on;
xlabel('X (km)');
ylabel('Y (km)');
title('卫星轨道 (近交点坐标系)');
axis equal;
box on;

plot(0, 0, 'ro', 'MarkerSize', 8, 'LineWidth', 2);
text(0, 0, ' 地球中心', 'FontSize', 10, 'VerticalAlignment', 'bottom');

plot(r_pqw_recovered(1, 1), r_pqw_recovered(2, 1), 'go', 'MarkerSize', 8, 'LineWidth', 2);
text(r_pqw_recovered(1, 1), r_pqw_recovered(2, 1), ' 起点', 'FontSize', 9, 'VerticalAlignment', 'bottom');

view(3);
camlight('left');
lighting gouraud;
hold off;
