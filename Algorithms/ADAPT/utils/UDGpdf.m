function f = UDGpdf(x, alpha_1, alpha_2, lambda_1, lambda_2)
% The PDF of the abs. value of the difference of two gamma
% distributions:
% Z = |X - Y| ~ UDG(alpha_1, alpha_2, lambda_1, lambda_2) where
% X ~ Gamma(alpha_1, lambda_1), Y ~ Gamma(alpha_2, lambda_2)

% Author: Pasquale Antonante
% email: antonap@mit.edu
% Date: 2021-01-06

fx = @(z) gampdf(z, alpha_1, lambda_1);
fy = @(z) gampdf(z, alpha_2, lambda_2);
pdf = @(t) integral(@(y) (fx(t+y) + fx(-t+y)).*fy(y), 0, Inf);
f = arrayfun(@(t) pdf(t), x);
f = f(:)';
end