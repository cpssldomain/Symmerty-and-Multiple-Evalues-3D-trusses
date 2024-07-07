function [y, dydx] = sym_fun_gen_cluster3(x,c)
y = x(1)*x(1)*x(2) + c*sin(x(2) + x(3));

dydx = [2*x(1)*x(2), x(1)*x(1)+c*cos(x(2)+x(3)), c*cos(x(2)+x(3))];
         
end