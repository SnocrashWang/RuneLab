function [para, omega_arr] = GN_alter(angle, time, omega)

% 默认omega初始化参数
if(~exist('omega','var'))
    omega = 1.942;  % 如果未出现该变量，则对其进行赋值
end

%% Gauss-Newton分步优化
    omega_arr = [omega];
    para = OLS(angle, time, omega);
    dw = 1e2;

    while abs(dw) > 1e-4
        syms w
        func = -para(1) / w * cos(w * time + para(3)) + (2.090 - para(1)) * time + para(4) - angle;
        Jac = jacobian(func, w);
        f = eval(subs(func, w, omega));
        J = eval(subs(Jac, w, omega));

        dw = -inv(J' * J) * J' * f';
        omega = omega + 3 * dw;
        omega_arr = [omega_arr, omega];
        para = OLS(angle, time, omega);
    end

end
