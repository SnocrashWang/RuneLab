function [time, ori, noisy, param] = getData(num, mode)

%% 判断数据模式
% 默认大符
if(~exist('mode','var'))
    mode = 0;  % 如果未出现该变量，则对其进行赋值
end

% 噪声指数
pose_noise = 0.;
gauss_noise = 0.05;
sparse_noise = [0.00, 0.5];    % density, sigma                                                 


% 加入噪声和稳态误差
rand_pose = randn(1,2) * pose_noise;

if mode == 0
    %% 大符
    % 旋转方程参数
    a = 0.780 + rand * (1.045 - 0.780);
    omega = 1.884 + rand * (2.000 - 1.884);
    % if omega_is_precise
    %     omega = 1.942;
    % end
    b = 2.090 - a;
    phi = rand * (2 * pi);

    param = [a, omega, phi, rand_pose];

    % 生成角度数据
    time = (0:num-1) / 100;
    angle_ori = -a / omega * cos(omega * time + phi) + b * time;

else
    %% 小符
    param = [10 * 2 * pi / 60, 0];
    % 生成角度数据
    time = (0:num-1) / 100;
    angle_ori = param(1) * time;
end

%% 后处理
% 加入噪声
angle_noisy = atan2(sin(angle_ori) * cos(rand_pose(1)), cos(angle_ori) * cos(rand_pose(2)));
angle_noisy = angle_noisy + randn(1, num) * gauss_noise;
angle_noisy = angle_noisy + sprandn(1, num, sparse_noise(1)) * sparse_noise(2);

% 归一化
angle_ori = angle_ori - angle_noisy(1);
angle_noisy = angle_noisy - angle_noisy(1);

for i = 2:num
    while angle_noisy(i) - angle_noisy(i-1) > pi
        angle_noisy(i) = angle_noisy(i) - 2 * pi;
    end
    while angle_noisy(i) - angle_noisy(i-1) < -pi
        angle_noisy(i) = angle_noisy(i) + 2 * pi;
    end
end

ori = angle_ori;
noisy = angle_noisy;


end
