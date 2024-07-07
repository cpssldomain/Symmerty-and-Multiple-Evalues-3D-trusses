function [xmax, dxmax_dx] = p_norm_stable(x, p)
x0 = max(x);
y = x / x0;
ypsum = sum(y.^p);
ymax = ypsum .^ (1.0 / p);
xmax = ymax * x0;
dxmax_dx = ypsum .^ (1.0 / p - 1.0) * (y .^ (p - 1.0));
         
end