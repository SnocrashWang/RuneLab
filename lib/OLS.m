function [x] = OLS(angle, time, omega)
%% 最小二乘法
    A = [sin(omega * time'), cos(omega * time'), time', ones(length(time), 1)];
    x = (A' * A) \ A' * angle';
end

