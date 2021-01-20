function [inliers, info] = greedy(problem, f, varargin)
%greedy - Greedy algorithm
% See: 
%    - "Outlier-Robust Estimation: Hardness, Minimally-Tuned Algorithms, and Applications" 
%      Antonante, Tzoumas, Yang, Carlone (2020).
%      https://arxiv.org/pdf/2007.15109.pdf
%
% Syntax:  [inliers, info] = greedy(problem, f, ...)
%
% Inputs:
%    problem - the stracture representing a generic problem (see EmptyProblem)
%    f - a function_handler of a non-minimal solver for problem
%
% Options:
%    - Formulation: Maximum-Consensus (MC) or Minimally-Trimmed-Square (MTS, default)
%    - CostThreshold: function_handle or scalar, if cost below the value, 
%         the alg. can stop (default: Inf)
%    - MinSolverMeasurements: minimum number of points that the non-minimal
%         solvers require to solve the problem (default: 1)
%
% Outputs:
%    inliers - the indices of the inliers
%    info - structure containing extended information about the xecution
% 
% Example:
%   problem = linearRegressionProblem(100, 0.8);
%   [inliers, info] = greedy(problem, @leastSquareNorm2, ...
%     'CostThreshold', @(n) chi2inv(0.99, n*problem.dof)*problem.MeasurementNoiseStd^2, ...
%     'MinSolverMeasurements', problem.MinMeasurements);

% Author: Pasquale Antonante
% email: antonap@mit.edu
% Date: 2021-01-06

assert(isProblem(problem), 'The problem doesn''t contain required fields');
start_t = now;
[inliers, info] = causalGreedy(problem, f, varargin{:});
end_t = now;
info.time = (end_t - start_t) * 24 * 60 * 60; % serial date number to sec
info.Algorithm = 'Greedy';
end

function [inliers, info] = causalGreedy(problem, f, varargin)
params = inputParser;
params.CaseSensitive = false;
params.addParameter('Formulation', 'MTS', @(x) ischar(x) && any(strcmpi({'MC', 'MTS'}, x)));
params.addParameter('CostThreshold', 1e-3, @(x) (isnumeric(x) && x>=0) || isa(x, 'function_handle'));
params.addParameter('MinSolverMeasurements', 1, @(x) isnumeric(x) && isscalar(x) && (x==floor(x)) && isfinite(x));
params.addParameter('Debug', false, @(x) islogical(x));
params.parse(varargin{:});

inliers = 1:problem.N;
try
    [f_val, f_info] = f(problem, 'Inliers', inliers);
catch err
    fprintf('Error message: %s', err.message)
    error("Could not run the global solver")
end
assert(isfield(f_info, 'residuals'), 'f should compute residuals');
switch params.Results.Formulation
    case 'MTS'
        cost = f_val;
    case 'MC'
        cost = max(f_info.residuals(inliers));
end
if isa(params.Results.CostThreshold, 'function_handle')
    cost_th = params.Results.CostThreshold(length(inliers));
else
    cost_th = params.Results.CostThreshold;
end

if params.Results.Debug
    S = {};
    S{1} = 1:problem.N;
    residual_hist = [];
    residual_hist(1,:) = f_info.residuals;
end

t = 1;
while cost > cost_th
    Z = inliers;
    T = setdiff(Z, problem.priors);
    [~, j] = max(f_info.residuals(T));
    t = t+1;
    inliers =  setdiff(Z, T(j));   
    try
        [f_val, f_info] = f(problem, 'Inliers', inliers);
    catch err
        fprintf('Error message: %s', err.message)
        error("Could not run the global solver")
    end
    switch params.Results.Formulation
        case 'MTS'
            cost = f_val;
        case 'MC'
            cost = max(f_info.residuals(inliers));
    end
    if isa(params.Results.CostThreshold, 'function_handle')
        cost_th = params.Results.CostThreshold(length(inliers));
    else
        cost_th = params.Results.CostThreshold;
    end
    if params.Results.Debug
        S{t} = inliers;
        residual_hist(t,:) = f_info.residuals;
    end
    if length(inliers) == params.Results.MinSolverMeasurements
        break
    end
end

info = struct();
info.Iterations = t;
info.stopping = 'Noise Threshold';
if params.Results.Debug
    info.SetHistory = S;
    info.ResidualsHistory = residual_hist;
end
end

