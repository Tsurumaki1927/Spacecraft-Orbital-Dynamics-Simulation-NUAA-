clear; clc; close all;

mu = 398600.4418; % 地球引力参数 (km^3/s^2) 
% 轨道参数（固定值） 
h = 58930;          % km^2/s
i = deg2rad(39.67); % 倾角
Omega = deg2rad(130.32); % RAAN
omega = deg2rad(42.373);  % 近地点幅角
e_values = [0.3, 0.4, 0.5, 0.6, 0.7, 0.8]; % 偏心率
colors = lines(length(e_values));

% 生成轨道轨迹点 
num_points = 200; %200points
ta = linspace(0, 2*pi, num_points); 

% 初始化图像
figure(1); set(gcf, 'Name', '2D'); 
figure(2); set(gcf, 'Name', '3D'); 

for k = 1:length(e_values) 
    e = e_values(k); 
    
    % 计算200个点
    r_all = zeros(num_points, 3); 
    for j = 1:num_points
        coe = [h, e, Omega, i, omega, ta(j)]; 
        [r, ~] = sv_from_coe(coe); 
        r_all(j, :) = r; 
    end
    
    % 二维
    figure(1); 
    subplot(2, 3, k); 
    plot(r_all(:,1), r_all(:,2), 'b', 'LineWidth', 1.5); 
    title(['e = ', num2str(e, '%.2f')]); 
    xlabel('X (km)'); ylabel('Y (km)'); 
    axis equal; grid on; 
    
    % 三维轨迹 
    figure(2); 
    subplot(2, 3, k); 
    hold on; 
    
    
    % 绘制轨道
    plot3(r_all(:,1), r_all(:,2), r_all(:,3), 'r', 'LineWidth', 1.5); 
    
    title(['3D e = ', num2str(e, '%.2f')]); 
    xlabel('X'); ylabel('Y'); zlabel('Z'); 
    axis equal; grid on; view(3);
    hold off; 
end



figure(3);
set(gcf, 'Name', '2D所有轨道');
hold on;

for k = 1:length(e_values)
    e = e_values(k); 
    
    % 计算轨道点
    r_all = zeros(num_points, 3); 
    for j = 1:num_points
        coe = [h, e, Omega, i, omega, ta(j)]; 
        [r, ~] = sv_from_coe(coe); 
        r_all(j, :) = r; 
    end
    
    % 绘制二维轨道 (X-Y平面)
    plot(r_all(:,1), r_all(:,2), 'Color', colors(k,:), 'LineWidth', 1.5);
end

% 添加参考线和标注
plot([0, 0], [-10000, 10000], 'k--', 'LineWidth', 0.5); % Y轴
plot([-10000, 10000], [0, 0], 'k--', 'LineWidth', 0.5); % X轴
xlabel('X (km)');
ylabel('Y (km)');
title('二维的轨道在同一图内的对比');
axis equal;
grid on;
legend(arrayfun(@(e) sprintf('e = %.1f', e), e_values, 'UniformOutput', false), ...
       'Location', 'bestoutside', 'FontSize', 8);
xlim([-120000, 120000]); % 自动调整范围
ylim([-120000, 120000]);

% ======== 三维轨道图 (所有e值在同一图) ========
figure(4);
set(gcf, 'Name', '3D所有轨道');
hold on;

for k = 1:length(e_values)
    e = e_values(k); 
    
    % 计算轨道点
    r_all = zeros(num_points, 3); 
    for j = 1:num_points
        coe = [h, e, Omega, i, omega, ta(j)]; 
        [r, ~] = sv_from_coe(coe); 
        r_all(j, :) = r; 
    end
    
    % 绘制三维轨道
    plot3(r_all(:,1), r_all(:,2), r_all(:,3), 'Color', colors(k,:), 'LineWidth', 1.5);
end

% 添加参考线和标注
plot3([0, 0], [0, 0], [-15000, 15000], 'k--', 'LineWidth', 0.5); % Z轴
plot3([0, 0], [-15000, 15000], [0, 0], 'k--', 'LineWidth', 0.5); % Y轴
plot3([-15000, 15000], [0, 0], [0, 0], 'k--', 'LineWidth', 0.5); % X轴
xlabel('X (km)');
ylabel('Y (km)');
zlabel('Z (km)');
title('三维的轨道在同一图内的对比');
axis equal;
grid on;
view(3);
legend(arrayfun(@(e) sprintf('e = %.1f', e), e_values, 'UniformOutput', false), ...
       'Location', 'bestoutside', 'FontSize', 8);
xlim([-12000, 12000]);
ylim([-12000, 12000]);
zlim([-12000, 12000]);





%%%%%%%%%%%%%%%%%%%%%
function [r, v] = sv_from_coe(coe) 

    mu = 398600.4418; 
    
    h    = coe(1); 
    e    = coe(2); 
    RA   = coe(3); 
    incl = coe(4); 
    w    = coe(5); 
    TA   = coe(6); 
    

    rp = (h^2 / mu) / (1 + e * cos(TA)) * [cos(TA); sin(TA); 0]; 
    vp = (mu / h) * [-sin(TA); e + cos(TA); 0]; 
    

    
    R3_RA = [ cos(RA),  sin(RA), 0; 
             -sin(RA),  cos(RA), 0; 
                   0,        0, 1]; 
                   
    R1_i = [1,        0,         0; 
            0,  cos(incl), sin(incl); 
            0, -sin(incl), cos(incl)]; 
            
    R3_w = [ cos(w),  sin(w), 0; 
            -sin(w),  cos(w), 0; 
                  0,       0, 1]; 
    

    Q_pX = R3_w * R1_i * R3_RA; 
    

    r = Q_pX' * rp; 
    v = Q_pX' * vp; 
    

    r = r'; 
    v = v'; 
end