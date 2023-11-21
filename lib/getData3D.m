function [time, ori, noisy, param] = getData3D(num, mode)

%% 判断数据模式
% 默认大符
if(~exist('mode','var'))
    mode = 0;  % 如果未出现该变量，则对其进行赋值
end

% 噪声指数
gauss_noise = 0.01;
sparse_noise = [0.00, 0.5];    % density, sigma

% 空间位置
theta = -pi + rand * (2 * pi);
dist_h = 6.5 + rand;
dist_v = 1.3 + rand * 2;
center = [dist_h * cos(theta) + (rand-0.5); dist_h * sin(theta) + (rand-0.5); dist_v];

if mode == 0
    %% 大符
    % 旋转方程参数
    a = 0.780 + rand * (1.045 - 0.780);
    omega = 1.884 + rand * (2.000 - 1.884);
    b = 2.090 - a;
    phi = -pi + rand * (2 * pi);
    c = -pi + rand * (2 * pi);

    param = [center; theta; a; omega; phi; c];

    % 生成三维数据
    time = (0:num-1) / 100;
    angle = -a / omega * cos(omega * time + phi) + b * time + c;

else
    %% 小符
    param = [center', theta, 10 * 2 * pi / 60, 0];
    % 生成三维数据
    time = (0:num-1) / 100;
    angle = param(1) * time;
end
r = 0.7;
position_ori = center + [r * cos(angle) * cos(theta + pi/2); r * cos(angle) * sin(theta + pi/2); r * sin(angle)];

%% 加入噪声
position_noisy = position_ori + randn(3, num) * gauss_noise;
position_noisy = position_noisy + sprandn(3, num, sparse_noise(1)) * sparse_noise(2);

ori = position_ori;
noisy = position_noisy;


end
