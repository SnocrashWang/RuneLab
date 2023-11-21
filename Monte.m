clc;
clear;
close all;
addpath(genpath('.\lib'));

%% 全局参数
dt = 0.01;
data_len = 600;
predict_len = 30;
monte = 10;
down_rate = 5;
% 0 for GN, 1 for GN_alter, 2 for FFT
method = 1

%% 数据记录
Para = zeros(4,monte);
Para_bar = zeros(4,monte);
Omega = zeros(1,monte);
Omega_bar = zeros(1,monte);
Dist = zeros(1,monte);
omega_arr = [0];
wrong_cnt = 0;

for i = 1:monte
    [time, angle_ori, angle_noisy, param] = getData(data_len+predict_len);
    spd_ori = (angle_ori(2:end) - angle_ori(1:end-1)) / dt;
    spd_noisy = (angle_noisy(2:end) - angle_noisy(1:end-1)) / dt;

    angle_input = downsample(meanFilter(angle_noisy(1:data_len), down_rate), down_rate);
    time_input = downsample(time(1:data_len), down_rate);
    spd_input = meanFilter(spd_noisy(1:data_len-1), down_rate);

    %% 拟合
    if method == 0              % 高斯牛顿法
        [para, omega_arr] = GN(angle_input, time_input);
        omega = omega_arr(end);
        % 迭代方向错误
        if (param(2)-omega)/(param(2)-omega_arr(1)) > 1
            wrong_cnt=wrong_cnt+1;
            fprintf("Wrong direction in optimizing\n")
%             figure;
%             plot((1:length(omega_arr)), omega_arr, "blue");
%             hold on;
%             plot((1:length(omega_arr)), ones(1, length(omega_arr)) * param(2), "green");
        end
    elseif method == 1          % 交替迭代高斯牛顿法
        [para, omega_arr] = GN_alter(angle_input, time_input);
        omega = omega_arr(end);
        % 迭代方向错误
        if (param(2)-omega)/(param(2)-omega_arr(1)) > 1
            wrong_cnt=wrong_cnt+1;
            fprintf("Wrong direction in optimizing\n")
%             figure;
%             plot((1:length(omega_arr)), omega_arr, "blue");
%             hold on;
%             plot((1:length(omega_arr)), ones(1, length(omega_arr)) * param(2), "green");
        end
    elseif method == 2          % 带窗傅里叶变换法
        omega = FFT(spd_input);
        para = OLS(angle_input, time_input, omega);
    end

    Para_bar(:,i) = para;
    Omega(i) = param(2);
    Omega_bar(i) = omega;
    Dist(i) = (-para(1) / para(2) * cos(para(2) * time(end) + para(3)) + (2.090 - para(1)) * time(end) + para(4) - angle_ori(end)) * 700;
    fprintf("In monte %3d, iteration: %3d, omega: %.4f -> %.4f, predict point err: %7.4f\n", i, length(omega_arr)-1, Omega(i), Omega_bar(i), Dist(i));
end

%% 数据统计
fprintf("RMSE of omega: %.6f, RMSE of predict point err: %.6f\n", RMSE(Omega - Omega_bar), RMSE(Dist));
if method == 0 || method == 1
    fprintf("%.2f%% of tests have a wrong optimizing direction\n", wrong_cnt/monte*100);
end

figure;
subplot(2,1,1);
h_omega = histogram(Omega-Omega_bar);
title("Histgram of omega error")
xlabel("Error(rad)")
ylabel("Count")

subplot(2,1,2);
h_angle = histogram(Dist);
title("Histgram of predict point error")
xlabel("Error(mm)")
ylabel("Count")
