function [x, omega_arr] = GN(angle, time, omega)

% 默认omega初始化参数
if(~exist('omega','var'))
    omega = 1.942;  % 如果未出现该变量，则对其进行赋值
end

%% Gauss-Newton分步优化
    omega_arr = [omega];
    x = OLS(angle, time, omega);
    dw = 1e2;

    while abs(dw) > 1e-4
        syms w
        func = x(1) * sin(w * time) + x(2) * cos(w * time) + x(3) * time + x(4) - angle;
        Jac = jacobian(func, w);
        f = eval(subs(func, w, omega));
        J = eval(subs(Jac, w, omega));

        dw = -inv(J' * J) * J' * f';
        omega = omega + 3 * dw;
        omega_arr = [omega_arr, omega];
        x = OLS(angle, time, omega);
    end

%     figure;
%     scatter(time, angle, 1, "red");
%     hold on;
%     angle_fit = x(1) * sin(omega * time) + x(2) * cos(omega * time) + x(3) * time + x(4);
%     plot(time, angle_fit, "blue");
%     pause;


end
