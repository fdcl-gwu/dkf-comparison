function y  =sinx_over_x(x)
eps = 1e-6;
if abs(x) < eps
    y =- x^10/39916800 + x^8/362880 - x^6/5040 + x^4/120 - x^2/6 + 1;
else
    y = sin(x)/x;
end
end