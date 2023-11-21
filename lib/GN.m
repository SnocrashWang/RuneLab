function [para, omega_arr] = GN(angle, time, omega)

% 默认omega初始化参数
if(~exist('omega','var'))
    omega = 1.942;  % 如果未出现该变量，则对其进行赋值
end

%% Gauss-Newton分步优化
    omega_arr = [omega];
    para = OLS(angle, time, omega);            % a, omega, phi, c
    dw = 1e2;

    syms a w phi c
    func = -a / w * cos(w * time + phi) + (2.090 - a) * time + c - angle;
    Jac = jacobian(func, [a, w, phi, c]);
    while norm(dw) > 1e-3
        f = eval(subs(func, [a, w, phi, c], para));
        J = eval(subs(Jac, [a, w, phi, c], para));

        dw = -inv(J' * J) * J' * f';
        para = para + 1 * dw';
        omega_arr = [omega_arr, para(2)];
    end

end
