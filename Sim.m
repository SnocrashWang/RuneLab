clc;
clear;
close all;
addpath(genpath('.\lib'));

%% 全局参数
dt = 0.01;
sim_time = 100;
data_len = 600;
predict_len = 30;
acc_rate = 10;
down_rate = 5;
% 0 for GNpure, 1 for GN, 2 for FFT
method = 0

%% 数据记录
Para = zeros(5,sim_time);
Para_bar = zeros(5,sim_time);
X_bar = zeros(4,sim_time);
Omega = zeros(1,sim_time);
Omega_bar = zeros(1,sim_time);
Dist = zeros(1,sim_time);
Angle = zeros(1,sim_time);
omega_arr = [0];
wrong_cnt = 0;
omega = 1.942;

%% 获取数据
[time, angle_ori, angle_noisy, param] = getData(sim_time*acc_rate+data_len+predict_len);
spd_ori = (angle_ori(2:end) - angle_ori(1:end-1)) / dt;
spd_noisy = (angle_noisy(2:end) - angle_noisy(1:end-1)) / dt;

for i = 1:sim_time
    cur_time = i*acc_rate;
    angle_input = downsample(meanFilter(angle_noisy(cur_time:cur_time+data_len), down_rate), down_rate);
    time_input = downsample(time(cur_time:cur_time+data_len), down_rate);
%     angle_input = angle_noisy(cur_time:cur_time+data_len);
%     time_input = time(cur_time:cur_time+data_len);
    spd_input = meanFilter(spd_noisy(cur_time:cur_time+data_len-1), down_rate);

    %% 拟合
    if method == 0              % 高斯牛顿法
        [x, omega_arr] = GNpure(angle_input, time_input, omega);
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
        [x, omega_arr] = GN(angle_input, time_input, omega);
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
        x = OLS(angle_input, time_input, omega);
    end

    X_bar(:,i) = x;
    Omega(i) = param(2);
    Omega_bar(i) = omega;
    pre_idx = cur_time+data_len+predict_len;
    Angle(i) = x(1)*sin(omega*time(pre_idx)) + x(2)*cos(omega*time(pre_idx)) + x(3)*time(pre_idx) + x(4);
    Dist(i) = (Angle(i) - angle_ori(pre_idx)) * 700;
    fprintf("At simulator time %4d ms, iteration: %3d, omega: %.4f -> %.4f, predict point err: %7.4f\n", cur_time+data_len, length(omega_arr)-1, Omega(i), Omega_bar(i), Dist(i));
end

%% 数据统计
fprintf("RMSE of omega: %.6f, RMSE of predict point err: %.6f\n", RMSE(Omega - Omega_bar), RMSE(Dist));
if method == 1
    AB = sqrt(X_bar(1,:).^2 + X_bar(2,:).^2) .* Omega_bar + X_bar(3,:);
    fprintf("RMSE of (a+b): %.6f\n", RMSE(AB - 2.090));
end
if method == 0 || method == 1
    fprintf("%.2f%% of tests have a wrong optimizing direction\n", wrong_cnt/sim_time*100);
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

figure;
plot(time(acc_rate+data_len+predict_len:sim_time*acc_rate+data_len+predict_len), ...
    angle_ori(acc_rate+data_len+predict_len:sim_time*acc_rate+data_len+predict_len), "green");
hold on;
plot(time(acc_rate+data_len+predict_len:sim_time*acc_rate+data_len+predict_len), ...
    angle_noisy(acc_rate+data_len+predict_len:sim_time*acc_rate+data_len+predict_len), "red");
hold on;
scatter(time(acc_rate+data_len+predict_len:acc_rate:sim_time*acc_rate+data_len+predict_len), Angle, 50, "blue", '+');
title("Rune function with time")
xlabel("time(dt)")
ylabel("Angle(rad)")

figure;
subplot(2,1,1);
scatter(time(acc_rate+data_len+predict_len:acc_rate:sim_time*acc_rate+data_len+predict_len), Omega_bar, 2, "blue");
hold on;
plot(time(acc_rate+data_len+predict_len:acc_rate:sim_time*acc_rate+data_len+predict_len), Omega, "green");
title("Omega error with time")
xlabel("time(dt)")
ylabel("Omega")

subplot(2,1,2);
scatter(time(acc_rate+data_len+predict_len:acc_rate:sim_time*acc_rate+data_len+predict_len), Dist, 2, "blue");
hold on;
plot(time(acc_rate+data_len+predict_len:acc_rate:sim_time*acc_rate+data_len+predict_len), zeros(1, sim_time), "black");
title("Predict point error with time")
xlabel("time(dt)")
ylabel("Error(mm)")
