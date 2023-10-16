function [x, omega_arr] = GNpure(angle, time, omega)

% 默认omega初始化参数
if(~exist('omega','var'))
    omega = 1.942;  % 如果未出现该变量，则对其进行赋值
end

%% Gauss-Newton分步优化
    omega_arr = [omega];
    x = OLS(angle, time, omega);
    para = [sqrt(x(1)^2 + x(2)^2)*omega;omega;atan2(x(2), x(1));x(4)];  % omega, a, phi, c
    dw = 1e2;

    syms a w phi c
    func = -a / w * cos(w * time + phi) + (2.090 - a) * time + c - angle;
    Jac = jacobian(func, [a, w, phi, c]);
    while norm(dw) > 1e-3
        f = eval(subs(func, [a, w, phi, c], [para(1), para(2), para(3), para(4)]));
        J = eval(subs(Jac, [a, w, phi, c], [para(1), para(2), para(3), para(4)]));

        dw = -inv(J' * J) * J' * f';
        para = para + 1 * dw;
        omega_arr = [omega_arr, para(2)];
    end

    x = [para(1)/para(2)*sin(para(3)), -para(1)/para(2)*cos(para(3)), 2.090-para(1), para(4)];

%     figure;
%     scatter(time, angle, 1, "red");
%     hold on;
%     angle_fit = -para(1) / para(2) * cos(para(2) * time + para(3)) + (2.090 - para(1)) * time + para(4);
%     plot(time, angle_fit, "blue");
%     pause;


end
