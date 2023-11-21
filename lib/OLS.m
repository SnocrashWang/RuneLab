function [para] = OLS(angle, time, omega)
%% 最小二乘法
    A = [sin(omega * time'), cos(omega * time'), time', ones(length(time), 1)];
    x = (A' * A) \ A' * angle';
    para = [sqrt(x(1)^2 + x(2)^2)*omega, omega, atan2(x(1), -x(2)), x(4)];
end

