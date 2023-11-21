function [para, para_arr] = GN3D(position, time, param)

%% 初始设置
% 优化变量：center_x, center_y, center_z, theta, a, omega, phi, c
para = getInit(position, time);
para_arr = [para];

% 目标函数
syms x y z th a w phi c
unknowns = [x y z th a w phi c];
unknowns_3d = [x y z th];
unknowns_angle = [a w phi c];
func_angle = -a / w * cos(w * time + phi) + (2.090 - a) * time + c;
func = (x + 0.7 * cos(func_angle) * cos(th + pi/2) - position(1)).^2 + ...
       (y + 0.7 * cos(func_angle) * sin(th + pi/2) - position(2)).^2 + ...
       (z + 0.7 * sin(func_angle) - position(3)).^2;
func_3d = x
Jac = jacobian(func, unknowns);

%% 交替迭代优化
dw = 1e2;
while norm(dw) > 1e-3
    f = eval(subs(func, unknowns, para'));
    J = eval(subs(Jac, unknowns, para'));

    dw = -(J' * J) \ (J' * f');
    para = para + 1 * dw;
    para'
    para_arr = [para_arr, para];
end

% figure;
% scatter(time, angle, 1, "red");
% hold on;
% angle_fit = -para(1) / para(2) * cos(para(2) * time + para(3)) + (2.090 - para(1)) * time + para(4);
% plot(time, angle_fit, "blue");
% pause;

end


function para_init = getInit(position, time)
%% 计算优化变量初值
% center_x, center_y, center_z, theta, a, omega, phi, c

% theta取xy平面的斜率最小二乘拟合
xOy = position(1:2, :);
ab = (xOy * xOy') \ (xOy * ones(size(position,2),1));
theta = atan2(ab(2), ab(1));
rot = [cos(theta),-sin(theta), 0;
       sin(theta), cos(theta), 0;
       0, 0, 1];

face_position = rot' * position;
flat_position = face_position(2:3,:);

% 优化求解最优中心
% 初值取全体坐标均值
center = mean(flat_position, 2);
syms c_x c_y
func_circle = (flat_position(1,:)-c_x).^2 + (flat_position(2,:)-c_y).^2 - 0.7^2;
Jac_circle = jacobian(func_circle, [c_x c_y]);

dw = 1e2;
while norm(dw) > 1e-4
    f = eval(subs(func_circle, [c_x c_y], center'));
    J = eval(subs(Jac_circle, [c_x c_y], center'));
    dw = -(J' * J) \ (J' * f');
    center = center + 1 * dw;
end

% 还原三维坐标
face_center = eval(subs([mean(face_position(1,:)); c_x; c_y], [c_x c_y], [center(1) center(2)]));
center = rot * face_center;

% 制作连续的角度数据
flat_circle = flat_position - face_center(2:3);
angle = atan2(flat_circle(2,:), flat_circle(1,:));
for i = 2:size(angle, 2)
    while angle(i) - angle(i-1) > pi
        angle(i) = angle(i) - 2 * pi;
    end
    while angle(i) - angle(i-1) < -pi
        angle(i) = angle(i) + 2 * pi;
    end
end

% % 旋转部分参数采用定omega的最小二乘法求取
% omega = 1.942;
% angle_para = OLS(angle, time, omega);
% a = sqrt(angle_para(1)^2 + angle_para(2)^2) * omega;
% phi = atan2(angle_para(1), -angle_para(2));
% c = angle_para(4);

[angle_para, ~] = GN(angle, time);

% figure
% scatter(flat_position(1,:), flat_position(2,:))
% xlabel('x')
% ylabel('y')
% plot(time, angle)

% para_init = [center'; theta; a; omega; phi; c]
para_init = [center; theta; angle_para];
para_init'
end

