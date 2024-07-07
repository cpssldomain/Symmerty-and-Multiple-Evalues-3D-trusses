function [y, dydx] = unsym_fun_gen(x)
y = x(1)*x(1)*x(2)*x(3) + x(2)*x(2)*x(3)*x(1) + x(3)*x(3)*x(1);

dydx = [2*x(1)*x(2)*x(3) + x(2)*x(2)*x(3) + x(3)*x(3),
        x(1)*x(1)*x(3) + 2*x(2)*x(3)*x(1) + 0,
        x(1)*x(1)*x(2) + x(2)*x(2)*x(1) + 2*x(3)*x(1)];
         
end