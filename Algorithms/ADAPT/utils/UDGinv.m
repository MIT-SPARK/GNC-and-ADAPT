function [v, p] = UDGinv(p, alpha_1, alpha_2, lambda_1, lambda_2)
% The inverse CDF of the abs. value of the difference of two gamma
% distributions:
% Z = |X - Y| ~ UDG(alpha_1, alpha_2, lambda_1, lambda_2) where
% X ~ Gamma(alpha_1, lambda_1), Y ~ Gamma(alpha_2, lambda_2)

% Author: Pasquale Antonante
% email: antonap@mit.edu
% Date: 2021-01-06

method = 'numeric';
N = 1e4;
[F, x] = UDGcdf(alpha_1, alpha_2, lambda_1, lambda_2, 'Method', method, 'Samples', N);
v = quadratic_interp(p, x, F);
end

function v = quadratic_interp(p, x, F)
i = find(diff(sign(F-p)))+1;
if length(i)>1
    v = x(i(1));
    return
end
Y = [F(i-1)-p; F(i)-p; F(i+1)-p];
A = [x(i-1)^2 x(i-1) 1; x(i)^2 x(i) 1; x(i+1)^2 x(i+1) 1];
C = A\Y;
a = C(1); b = C(2); c = C(3);
v = (-b + sqrt(b^2 - 4*a*c)) / (2*a);
end

function v = linear_interp(p, x, F)
i = find(diff(sign(F-p)))+1;
v = (x(i)-x(i-1))/(F(i)-F(i-1))*(p-F(i-1))+x(i-1);
end
