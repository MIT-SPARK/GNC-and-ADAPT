function [inliers, info] = gnc(problem, f, varargin)
%GNC - Graduated Non-Convexity
% Implementation of: 
%    - "Graduated Non-Convexity for Robust Spatial Perception: From Non-Minimal Solvers to Global Outlier Rejection"
%      Yang, Antonante, Tzoumas, Carlone (2020). IEEE Robotics and Automation Letters (RA-L), 5(2), 1127–1134.
%      https://arxiv.org/pdf/1909.08605.pdf
%    - "Outlier-Robust Estimation: Hardness, Minimally-Tuned Algorithms, and Applications" 
%      Antonante, Tzoumas, Yang, Carlone (2020).
%      https://arxiv.org/pdf/2007.15109.pdf
%
% Syntax:  [inliers, info] = gnc(problem, f, ...)
%
% Inputs:
%    problem - the stracture representing a generic problem (see EmptyProblem)
%    f - a function_handler of a non-minimal solver for problem
%
% Options:
%    - ContinuationFactor: continuation factor for mu
%    - InitialMu: initial mu (default: auto)
%    - NoiseBound: inlier threshold
%    - MaxIterations: maximum number of iterations
%    - FixPriorsWeights: fix priors weights to 1 across all iterations
%    - CostThreshold: if weighted sum of squared residuals is below this value, GNC stops
%    - Debug: whether or not to enable the debug information
%
% Outputs:
%    inliers - the indices of the inliers
%    info - structure containing extended information about the xecution
% 
% Example:
%   problem = linearRegressionProblem(100, 0.8);
%   [inliers, info] = gnc(problem, @leastSquareNorm2, ...
%       'NoiseBound', chi2inv(0.99, problem.dof)*problem.MeasurementNoiseStd^2);

% Author: Pasquale Antonante
% email: antonap@mit.edu
% Date: 2021-01-06

params = inputParser;
params.KeepUnmatched = true;
params.addParameter('NoiseBound', 0, @(x) isscalar(x) && x>=0 && isfinite(x));
params.parse(varargin{:});
assert(isProblem(problem), 'The problem doesn''t contain required fields.');

start_t = now;
if params.Results.NoiseBound > 0
    [inliers, info] = gnc_vanilla(problem, f, varargin{:});
else
    error('You need to set the Noise Bound')
end
end_t = now;
info.time = (end_t - start_t) * 24 * 60 * 60; % serial date number to sec
end

function [inliers, info] = gnc_vanilla(problem, f, varargin)
params = inputParser;
params.KeepUnmatched = true;
params.addParameter('ContinuationFactor', 1.4, @(x) isscalar(x) && x>0);
params.addParameter('InitialMu', 'auto', @(x) ischar(x) || (isscalar(x) && x>0));
params.addParameter('NoiseBound', 0, @(x) isscalar(x) && x>=0 && isfinite(x));
params.addParameter('MaxIterations', 1e3, @(x) isscalar(x));
params.addParameter('FixPriorsWeights', true, @(x) islogical(x));
params.addParameter('CostThreshold', 0, @(x) isscalar(x));
params.addParameter('Debug', false, @(x) islogical(x));
params.addParameter('init_', 0);
params.parse(varargin{:});

max_iterations = params.Results.MaxIterations;

if params.Results.Debug
    residuals_history = [];
    weights_history = [];
end

if ismember('init_', params.UsingDefaults)
    try
        [~, f_info] = f(problem);
    catch err
        fprintf('Error message: %s', err.message)
        error("Could not run the global solver")
    end
    assert(isfield(f_info, 'residuals'), 'f should compute residuals');
else
    f_info = params.Results.('init_');
end
barc2 = params.Results.NoiseBound;
weights = ones(1, problem.N);
if ischar(params.Results.InitialMu) && strcmpi(params.Results.InitialMu, 'auto')
    mu = 1 / (2 * max(f_info.residuals) / barc2 - 1 );
elseif isnumeric(params.Results.InitialMu)
    mu = params.Results.InitialMu;
else
    error('Invalid value for InitialMu')
end
    
prev_f_cost = 0;

i = 1;
info.stopping = 'MaxIterations'; % Worst case stopping
info.status = 1; % worst case status: failure
while i < max_iterations
    if params.Results.Debug
        residuals_history(i,:) = f_info.residuals(:)';
        weights_history(i, :) = weights;
    end
    weights = gncWeightsUpdate(weights, mu, f_info.residuals, barc2);
    if params.Results.FixPriorsWeights
        weights(problem.priors) = 1;
    end
    try
        [~, f_info] = f(problem, 'Weights', weights);
    catch err
        fprintf('Error message: %s', err.message)
        error("Could not run the global solver")
    end
    f_cost = sum(f_info.residuals(:) .* weights(:));
    cost_diff = abs(f_cost - prev_f_cost);
    prev_f_cost = f_cost;
    
    if (cost_diff < params.Results.CostThreshold) || areBinaryWeights(weights) 
        if params.Results.Debug
            residuals_history(i+1,:) = f_info.residuals(:)';
            weights_history(i+1, :) = weights;
        end
        info.residuals = f_info.residuals(:)';
        info.stopping = 'CostThreshold';
        info.status = 1;
        break       
    end
    mu = mu * params.Results.ContinuationFactor;
    i = i + 1;
end
inliers = find(weights>1-eps);
info.Iterations = i;
info.params = params.Results;
info.Algorithm = 'GNC';
if params.Results.Debug
    info.mu = mu;
    info.barc2History = repmat(barc2, i+1, 1);
    info.ResidualsHistory = residuals_history;
    info.WeightsHistory = weights_history;
    info.CostDiff = cost_diff;
    info.AreBinaryWeights = areBinaryWeights(weights);
end
end