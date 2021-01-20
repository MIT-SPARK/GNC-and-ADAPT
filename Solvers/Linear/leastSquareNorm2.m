function [f_val, info] = leastSquareNorm2(problem, varargin)
% Solve a problem of the type $$\min_x \sum_{i\in\text{selections}} \|(A_i*x - y_i\|^2$$
% See also: linearRegressionProblem

assert(isProblem(problem) && strcmpi(problem.type, 'linear-estimation'), ...
    'You cannot solve this problem with linear least square')
measurements =  1:problem.N;
default_weights = ones(1, problem.N);

params = inputParser;
params.CaseSensitive = false;
params.addParameter('Inliers', measurements, @(x) isnumeric(x) && all(floor(x)==x));
params.addParameter('Weights', default_weights, @(x) isnumeric(x));
params.addParameter('Estimate', [], @(x) isnumeric(x));
params.parse(varargin{:});
assert(numel(params.Results.Weights) == problem.N, 'The weights should be a 1xN vector');

inliers = params.Results.Inliers(:)';
w = params.Results.Weights(:);

% Build weighted matrices
A_ = cell(1, length(inliers));
for i = inliers
    A_{i} = sqrt(w(i)).*problem.A{i};
end
A_ = sparse(vertcat(A_{:}));
y_ = cell(1, length(inliers));
for i = inliers
    y_{i} = sqrt(w(i)).*problem.y{i};
end
y_ = sparse(vertcat(y_{:}));

if isempty(params.Results.Estimate)
    x = (A_'*A_)\(A_'*y_);
else
    x = params.Results.Estimate;
end
f_val = norm(A_*x-y_)^2;

% Residues
residues = zeros(problem.MeasurementDimension, problem.N);
for i = 1:problem.N
    residues(:,i) = problem.A{i}*x(:) - problem.y{i};
end
% Residuals
residuals = zeros(1, problem.N);
for i = 1:problem.N
    residuals(i) = norm(residues(:,i)).^2;
end

info = struct('x', x, 'residuals', residuals, 'residues', residues');
end
