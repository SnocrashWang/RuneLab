clc;
clear;
close all;
addpath(genpath('.\lib'));

%% 获取数据集
num = 400;
[time, ori, noisy, param] = getData3D(num);
% center_x, center_y, center_z, theta, a, omega, phi, c
param'

% 绘图
figure
scatter3(ori(1,:), ori(2,:), ori(3,:))
hold on
scatter3(noisy(1,:), noisy(2,:), noisy(3,:))
% hold on
% scatter3(param(1), param(2), param(3))
xlabel('x')
ylabel('y')
zlabel('z')
% xlim([-8,8]);
% ylim([-8,8]);
% zlim([-8,8]);

% 降采样
down_rate = 5;
position_input(1,:) = downsample(meanFilter(noisy(1,:), down_rate), down_rate);
position_input(2,:) = downsample(meanFilter(noisy(2,:), down_rate), down_rate);
position_input(3,:)= downsample(meanFilter(noisy(3,:), down_rate), down_rate);
time_input = downsample(time, down_rate);

%% 求解运动方程
[x, para_arr] = GN3D(position_input, time_input, param);

