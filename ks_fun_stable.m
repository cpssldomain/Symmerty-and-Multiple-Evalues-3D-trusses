function [y, dydx] = ks_fun_stable(x, q)
x1 = max(x);
xexp_sum = sum(exp(q * (x-x1)));

y = x1 + 1 / q*log(xexp_sum);
dydx = 1.0 / xexp_sum * exp(q*(x-x1));
         
end