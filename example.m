%% SETTINGS - change here
num_samples = 100; % how many measurements are available
outliers_percentage = 0.8;

% To disable the execution of any algorithm feel free to comment out the
% relative section. `displayResults` will show only available results.

%% Problem generation
problem = linearRegressionProblem(num_samples, outliers_percentage);
problem.f = @leastSquareNorm2;
problem.MinMeasurements = 1;
printProblemSummary(problem);

% Thresholds
cfg.epsilon = chi2inv(0.99, problem.dof)*problem.MeasurementNoiseStd^2;
cfg.epsilon_mts = @(n) chi2inv(0.99, n*problem.dof)*problem.MeasurementNoiseStd^2;
cfg.epsilon_dmts = @(n1, n2) UDGinv(0.99, ...
    n1*problem.dof/2, n2*problem.dof/2, ... % shape parameters
    2*problem.MeasurementNoiseStd^2, 2*problem.MeasurementNoiseStd^2); % scale parameters

% Solvers
results = struct();
fprintf('Solving with:\n')

%% Solve -- Non-Robust (uses all measurements, including outliers)
fprintf('  - Least-Square\n')
start_t = now;
f_val = leastSquareNorm2(problem, 'Inliers', 1:problem.N);
runtime = (now-start_t)*24*60*60;

results.ls.algname = 'Least-Square';
results.ls.f_val = f_val;
[results.ls.stats.tp, results.ls.stats.fp] = detectionStats(problem, 1:problem.N);
results.ls.iterations = 1;
results.ls.time = runtime;

%% Solve -- Grount-Truth (only uses ground truth inliers)
fprintf('  - Ground-Truth Least-Square...\n')
inliers_gt = setdiff(1:problem.N, problem.outliers);
start_t = now;
f_val = leastSquareNorm2(problem, 'Inliers', inliers_gt);
runtime = (now-start_t)*24*60*60;

results.gt.algname = 'Ground-Truth';
results.gt.f_val = f_val;
[results.gt.stats.tp, results.gt.stats.fp] = detectionStats(problem, inliers_gt);
results.gt.iterations = 1;
results.gt.time = runtime;

%% Solve -- RANSAC
fprintf('  - RANSAC\n')
start_t = now;
[model, ransac_inliers, ransac_itr] = ransac(problem.N, ...
    problem, @ransacFitFcn, @ransacResidualsFcn, ...
    problem.MinMeasurements, cfg.epsilon, ...
    1.0, 1000);
runtime = (now-start_t)*24*60*60;

results.ransac.algname = 'RANSAC';
results.ransac.f_val = leastSquareNorm2(problem, 'Inliers', ransac_inliers);
[results.ransac.stats.tp, results.ransac.stats.fp] = ...
    detectionStats(problem, ransac_inliers);
results.ransac.iterations = ransac_itr;
results.ransac.time = runtime;

%% Solve -- Greedy (Minimally-Trimmed Square)
fprintf('  - GREEDY (MTS)\n')
[inliers, info] = greedy(problem, problem.f, ...
    'Formulation', 'MTS', ...
    'CostThreshold', cfg.epsilon_mts, ...
    'MinSolverMeasurements', problem.MinMeasurements);

results.greedy_mts.algname = 'Greedy (MTS)';
results.greedy_mts.f_val = leastSquareNorm2(problem, 'Inliers', inliers);
[results.greedy_mts.stats.tp, results.greedy_mts.stats.fp] = ...
    detectionStats(problem, inliers);
results.greedy_mts.iterations = info.Iterations;
results.greedy_mts.time = info.time;

%% Solve -- ADAPT (MTS)
fprintf('  - ADAPT\n')
[inliers, info] = adapt(problem, problem.f, ...
    'CostThreshold',  cfg.epsilon_mts, ...
    'ConvergenceThreshold',  cfg.epsilon_dmts, ...
    'MinSolverMeasurements', problem.MinMeasurements);

results.adapt_mts.algname = 'ADAPT (MTS)';
results.adapt_mts.f_val = leastSquareNorm2(problem, 'Inliers', inliers);
[results.adapt_mts.stats.tp, results.adapt_mts.stats.fp] = ...
    detectionStats(problem, inliers);
results.adapt_mts.iterations = info.Iterations;
results.adapt_mts.time = info.time;

%% Solve -- GNC
fprintf('  - GNC\n')
[inliers, info] = gnc(problem, problem.f, 'NoiseBound', cfg.epsilon);

results.gnc.algname = 'GNC';
results.gnc.f_val = leastSquareNorm2(problem, 'Inliers', inliers);
[results.gnc.stats.tp, results.gnc.stats.fp] = ...
    detectionStats(problem, inliers);
results.gnc.iterations = info.Iterations;
results.gnc.time = info.time;

%% Results
fprintf('Results:\n')
displayResults(results);

%% ------------------------------------
function displayResults(results)
err = []; tp = []; fp = []; time = []; itr = []; alg_name = {};
field_names = fieldnames(results);
for k=1:numel(field_names)
    err(end+1) = results.(field_names{k}).f_val;
    tp(end+1) = results.(field_names{k}).stats.tp;
    fp(end+1) = results.(field_names{k}).stats.fp;
    time(end+1) = results.(field_names{k}).time;
    itr(end+1) = results.(field_names{k}).iterations;
    alg_name{end+1} = results.(field_names{k}).algname;
end
T = table(err', 100.*tp', 100.*fp', itr', time', ...
    'RowNames', alg_name, ...
    'VariableNames', {'Error', 'TP', 'FP', 'Itr', 'Time'});
disp(T)
end

function model = ransacFitFcn(problem, selection)
[~, model] = leastSquareNorm2(problem, 'Inliers', selection);
end

function residuals = ransacResidualsFcn(~, model)
residuals = model.residuals;
end
