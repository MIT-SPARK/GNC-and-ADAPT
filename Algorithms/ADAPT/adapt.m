function [inliers, info] = adapt(problem, f, varargin)
%ADAPT - Adaptive Trimming algorithm
% Implementation of: 
%    - "Outlier-Robust Spatial Perception: Hardness, General-Purpose Algorithms, and Guarantees."
%      Tzoumas, Antonante, Carlone (2019). IEEE/RSJ Intl. Conf. on Intelligent Robots and Systems (IROS).
%      https://arxiv.org/pdf/1903.11683.pdf
%    - "Outlier-Robust Estimation: Hardness, Minimally-Tuned Algorithms, and Applications" 
%      Antonante, Tzoumas, Yang, Carlone (2020).
%      https://arxiv.org/pdf/2007.15109.pdf
%
% Syntax:  [inliers, info] = adapt(problem, f, ...)
%
% Inputs:
%    problem - the stracture representing a generic problem (see EmptyProblem)
%    f - a function_handler of a non-minimal solver for problem
%
% Options:
%    - Formulation: Maximum-Consensus (MC) or Minimally-Trimmed-Square (MTS, default)
%    - CostThreshold: function_handle or scalar, if cost below the value, 
%         the alg. can stop (default: Inf)
%    - ConvergenceThreshold: function handle or scalar cost threshold, if
%         cost below the value, the alg. increments the convergence counter
%         (default: Inf)
%    - ConvergenceIterations: number of iterations before convergence
%         (default: 3)
%    - DiscountFactor: discount factor for the inlier threshold 
%         (default: 0.99)
%    - MaxIterations: maximum number of iterations
%    - MinSolverMeasurements: minimum number of points that the non-minimal
%         solvers require to solve the problem (default: 1)
%    - SelfCorrection: whether or not to enable the self correction mechanims
%    - Debug: whether or not to enable the debug information
%
% Outputs:
%    inliers - the indices of the inliers
%    info - structure containing extended information about the xecution
% 
% Example:
%    problem = linearRegressionProblem(100, 0.8);
%    [inliers, info] = adapt(problem, @leastSquareNorm2, ...
%        'CostThreshold',  chi2inv(0.99, problem.dof)*problem.MeasurementNoiseStd^2);

% Author: Pasquale Antonante
% email: antonap@mit.edu
% Date: 2021-01-06

assert(isProblem(problem), 'The problem doesn''t contain required fields');
params = inputParser;
params.CaseSensitive = false;
params.KeepUnmatched = true;
params.parse(varargin{:});

start_t = now;
[inliers, info] = adaptiveTrimming(problem, f,  varargin{:});
end_t = now;
info.time = (end_t - start_t) * 24 * 60 * 60; % serial date number to sec
end

%% Adapptive Trimming
function [inliers, info] = adaptiveTrimming(problem, f, varargin)
params = inputParser;
params.CaseSensitive = false;
params.KeepUnmatched = true;
% Vanilla (can be used in combination with autotuning)
params.addParameter('Formulation', 'MTS', @(x) ischar(x) && any(strcmpi({'MC', 'MTS'}, x)));
params.addParameter('CostThreshold', 1e-3, @(x) (isnumeric(x) && x>=0) || isa(x, 'function_handle'));
params.addParameter('ConvergenceThreshold', Inf, @(x) (isnumeric(x) && x>=0) || isa(x, 'function_handle'));
params.addParameter('ConvergenceIterations', 3, @(x) isnumeric(x) && isscalar(x) && (x==floor(x)) && isfinite(x) && x>=1);
params.addParameter('DiscountFactor', 0.99, @(x) isnumeric(x) && isscalar(x) && x>0 && x<=1);
params.addParameter('MaxIterations', 1e4, @(x) isnumeric(x) && x>0 && x==floor(x) && isfinite(x));
params.addParameter('MinSolverMeasurements', 1, @(x) isnumeric(x) && isscalar(x) && (x==floor(x)) && isfinite(x));
params.addParameter('SelfCorrection', true, @(x) islogical(x));
params.addParameter('Debug', false, @(x) islogical(x));
params.parse(varargin{:});

% Initialization
alpha = params.Results.DiscountFactor;
inliers = 1:problem.N;
try
    [cost, f_info] = f(problem, 'Inliers', inliers);
catch err
    fprintf('Error message: %s', err.message)
    error("Could not run the global solver")
end
assert(isfield(f_info, 'residuals'), 'f should compute residuals');
barc2 = alpha*max(f_info.residuals);
S{1} = inliers;
if params.Results.Debug
    barc2_history = barc2;
    residuals_history = f_info.residuals;
    cost_history(1) = -Inf;
end
last_cost_mts = cost;
last_cost_mc = max(f_info.residuals);

% Main loop
i = 1;
convergence_counter = 0;
stop = 'Max Iterations';
while i < params.Results.MaxIterations
    i = i+1;
    % Find outliers set
    if params.Results.SelfCorrection
        inliers = union(find(f_info.residuals < barc2), problem.priors);
    else
        Z = setdiff(inliers, problem.priors);
        [~, j] = max(f_info.residuals(Z));
        inliers = union(setdiff(inliers, Z(j)), problem.priors);
    end
    S{i} = inliers;
    % Update estimates
    try
        [f_mts, f_info] = f(problem, 'Inliers', inliers);
    catch err
        fprintf('Error message: %s', err.message)
        error("Could not run the global solver")
    end
    f_mc = max(f_info.residuals(inliers));
    % Costs computation
    switch params.Results.Formulation
        case 'MC'
            cost = f_mc;
            conv_cost = abs(last_cost_mc - f_mc);
        case 'MTS'
            cost = f_mts;
            conv_cost = abs(last_cost_mts - f_mts);
    end
    last_cost_mts = f_mts;
    last_cost_mc = f_mc;
    % Thresholds computation
    if isa(params.Results.CostThreshold, 'function_handle')
        switch params.Results.Formulation
            case 'MC'
                cost_th = params.Results.CostThreshold(); 
            case 'MTS'
                cost_th = params.Results.CostThreshold(length(inliers));
        end
    else
        cost_th = params.Results.CostThreshold;
    end
    % Debug info
    if params.Results.Debug
        barc2_history(i) = barc2;
        residuals_history(i,:) = f_info.residuals;
        cost_history(i) = conv_cost;
    end
    % Stopping criteria
    if cost < cost_th
        % Compute convergence threshold
        if isa(params.Results.ConvergenceThreshold, 'function_handle')
            n_curr = length(S{i});
            n_prev = length(S{i-1});
            conv_th = params.Results.ConvergenceThreshold(n_curr, n_prev); 
        else
            conv_th = params.Results.ConvergenceThreshold;
        end
        if conv_cost < conv_th
            convergence_counter = convergence_counter + 1;
            if convergence_counter == params.Results.ConvergenceIterations
                stop = 'Cost Threshold';
                inliers = S{i-params.Results.ConvergenceIterations+1};
                break
            end
        end
    else
        convergence_counter = 0;
    end
    % Update barc2
    barc2 = alpha*min([barc2, max(f_info.residuals(setdiff(inliers, problem.priors)))]);
    if length(inliers) < params.Results.MinSolverMeasurements
        info.stopping = 'MinMeasurements';
        break
    end
end
info.params = params.Results;
info.Iterations = i;
info.Algorithm = 'ADAPT';
info.stopping = stop;
if params.Results.Debug
    info.SetHistory = S;
    info.ResidualsHistory = residuals_history;
    info.barc2History = barc2_history;
    info.CostHistory = cost_history;
end
end

