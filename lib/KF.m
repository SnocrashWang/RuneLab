function [x, P] = KF(x, P, z, b, w)
%% 卡尔曼滤波
    % 模型
    dt = 0.01;
    if length(x) == 3
        F = [1, dt, 0;
             -w^2 * dt, 1, 0;
             0, 0, 1];
        G = [0; 0; b * dt];
        H = [1, 0, 1];
        Q = [1e-3, 0, 0;
             0, 1e-3, 0;
             0, 0, 5e-4];
        R = 1e-2;
    elseif length(x) == 2
        F = [1, dt;
             -w^2 * dt, 1];
        G = [0; 0];
        H = [1, 0];
        Q = [1e-3, 0;
             0, 1e-3];
        R = 1e-2;
    end

    % 推算
    x_ = F * x + G;
    P_ = F * P * F' + Q;
    K = P_ * H' / (H * P_ * H' + R);
    x = x_ + K * (z - H * x_);
    P = (eye(length(x)) - K * H) * P;
end

