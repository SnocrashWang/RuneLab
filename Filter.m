clc;
clear;
close all;
addpath(genpath('.\lib'));

%% 全局参数
data_len = 800;
ols_len = 300;
init_len = 200;
omega = 1.942;
state_num = 2

%% 数据
[time, angle_ori, angle_noisy, param_ori] = getData(data_len);
param_true = OLS(angle_ori, time, param_ori(2));
param_init = OLS(angle_noisy(1:init_len), time(1:init_len), omega);
angle_ref = [param_true(1)*sin(param_ori(2)*time); param_true(2)*cos(param_ori(2)*time); param_true(3)*time;param_true(4)*ones(1, data_len)];

if state_num == 3
    x_arr = zeros(3, data_len-init_len);
    x_arr(:,1) = [param_init(1)*sin(omega*time(init_len+1)) + param_init(2)*cos(omega*time(init_len+1));
                  omega*param_init(1)*cos(omega*time(init_len+1)) - omega*param_init(2)*sin(omega*time(init_len+1));
                  param_init(3)*time(init_len+1)+param_init(4)];
elseif state_num == 2
    x_arr = zeros(2, data_len-init_len);
    x_arr(:,1) = [param_init(1)*sin(omega*time(init_len+1)) + param_init(2)*cos(omega*time(init_len+1));
                  omega*param_init(1)*cos(omega*time(init_len+1)) - omega*param_init(2)*sin(omega*time(init_len+1))];
end
P = eye(state_num) * 1e-6;

for i = (init_len+2):data_len
    ols_idx = max(1,i-ols_len-1);
    param_init = OLS(angle_noisy(ols_idx:i), time(ols_idx:i), omega);
    if state_num == 3
        [x_arr(:,i-init_len), P] = KF(x_arr(:,i-init_len-1), P, angle_noisy(i), param_init(3), omega);
    elseif state_num == 2
        [x_arr(:,i-init_len), P] = KF(x_arr(:,i-init_len-1), P, angle_noisy(i)-param_init(3)*time(i)-param_init(4), param_init(3), omega);
    end
end

%% 数据统计
figure;
plot(time, angle_ori, "green");
hold on;
plot(time, angle_noisy, "red");

plot(time(init_len+1:data_len), x_arr(1,:)+param_true(3)*time(init_len+1:data_len)+param_true(4)*ones(1,data_len-init_len), "blue");
plot(time(init_len+1:data_len), x_arr(1,:), "m");
plot(time(init_len+1:data_len), x_arr(2,:), "y");
if state_num == 3
    plot(time(init_len+1:data_len), x_arr(3,:), "cyan");
end

plot(time, angle_ref(1,:)+angle_ref(2,:), "m-.");
plot(time, angle_ref(3,:)+angle_ref(4,:), "cyan--");

title("Rune function with time")
xlabel("time(dt)")
ylabel("Angle(rad)")
if state_num == 3
    legend("angle origin", "angle noizy", ...
           "angle KF", "KF state_1", "KF state_2", "KF state_3", ...
           "state_1 origin", "state_3 origin", ...
           "Location", "NorthWest")
elseif state_num == 2
    legend("angle origin", "angle noizy", ...
           "angle KF", "KF state_1", "KF state_2", ...
           "state_1 origin", "state_3 origin", ...
           "Location", "NorthWest")
end