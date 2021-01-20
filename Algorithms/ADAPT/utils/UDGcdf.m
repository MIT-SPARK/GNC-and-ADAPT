function [F,x] = UDGcdf(alpha_1, alpha_2, lambda_1, lambda_2, varargin)
% The CDF of the abs. value of the difference of two gamma
% distributions:
% Z = |X - Y| ~ UDG(alpha_1, alpha_2, lambda_1, lambda_2) where
% X ~ Gamma(alpha_1, lambda_1), Y ~ Gamma(alpha_2, lambda_2)

% Author: Pasquale Antonante
% email: antonap@mit.edu
% Date: 2021-01-12

params = inputParser;
params.CaseSensitive = false;
params.addParameter('Method', 'auto', @(x) ischar(x) && any(strcmpi({'Auto', 'Numeric', 'Integration'}, x)));
params.addParameter('Samples', [], @(x) isnumeric(x) && isscalar(x) && (x==floor(x)) && isfinite(x) && x > 1);
params.addParameter('X', []);
params.parse(varargin{:});

method = params.Results.Method;
if strcmpi('auto',method)
    if alpha_1 > 300 || alpha_2 > 300
        method = 'numeric';
    else
        method = 'integration';
    end
end

switch lower(method)
    case 'integration'
        tol = 5e-4;
        if isempty(params.Results.Samples)
            N = 100;
        else
            N = params.Results.Samples;
        end
        if isempty(params.Results.X)
            x_end = max([gaminv(1-tol, alpha_1, lambda_1), gaminv(1-tol, alpha_2, lambda_2)]);
            x = linspace(0, x_end, N);
        else
            x = params.Results.X;
        end
        f = UDGpdf(x, alpha_1, alpha_2, lambda_1, lambda_2);
        F = cumsim(x,f);
        F = F./F(end);
    case 'numeric'
        if isempty(params.Results.Samples)
            N = 5e4;
        else
            N = params.Results.Samples;
        end
        X = gamrnd(alpha_1, lambda_1, N, 1);
        Y = gamrnd(alpha_2, lambda_2, N, 1);
        Z = abs(X - Y); 
        [F, x] = ecdf(Z);
end

end
