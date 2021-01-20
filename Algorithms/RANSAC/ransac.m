function [model, inlier_idx, itr] = ransac(N, data, fit_fun, residuals_fun, ...
sample_size, inlier_threshold, inlier_ratio_threshold, max_iterations)
%ransac - RANSAC (RANdom SAmple Consensus) algorithm
%
% Syntax:  [output1,output2] = function_name(input1,input2,input3)
%
% Inputs:
%    N - Total nnumber of samples
%    data - The problem data (used by the fit_fcn and resudual_fcn)
%    fit_fun - A function that compute the model given the problem and a selection of samples
%    residuals_fun - A function that computes the resuduals given a problem and a model
%    sample_size - The number of samples required to compute the model
%    inlier_threshold - Threshold to consider a measurement an inlier
%    inlier_ratio_threshold - Expected inlier ratio (if available, otherwise 1)
%    max_iterations - Max number of iterations
%
% Outputs:
%    model - Description
%    inlier_idx - Description
%    itr - Description
%
% Function prototypes:
%    model = fit_fun(data, selection)
%    residuals = residuals_fun(data, model)
%
% Reference: Hartley, R., & Zisserman, A. (2003). Multiple view geometry in computer vision.
%
% Example:
% problem = linearRegressionProblem(100, 0.8);
% [model, inlier_idx, itr] = ransac(problem.N, ...
%     problem, ...
%     @ransacFitFcn, ...
%     @ransacResidualsFcn, ...
%     problem.MinMeasurements, ...
%     chi2inv(0.99, problem.dof)*problem.MeasurementNoiseStd^2, ...
%     1.0, 1000);
% 
% function model = ransacFitFcn(problem, selection)
% [~, model] = leastSquareNorm2(problem, 'Inliers', selection);
% end
% 
% function residuals = ransacResidualsFcn(problem, model)
% residuals = model.residuals;
% end
%

% Author: Pasquale Antonante
% email: antonap@mit.edu
% Date: 2021-01-12

assert(sample_size>0, 'The sample size must be greater than 0')
assert(N>0 && N>sample_size, 'N must be greater than 0 and than the sample size')
assert(max_iterations > 0, 'The maximum number of iteration must be greater than 0')
assert(inlier_ratio_threshold > 0 && inlier_ratio_threshold <= 1, ...
    'The inlier ratio threshold must be a number in (0,1)')

best_model = [];
best_inlier_ratio = -inf;
best_inliers = [];
adapt_max_itr = 1;
itr = 0;
while(itr < adapt_max_itr)
    rand_idx = datasample(1:N, sample_size, 'Replace',false);
    try
        model = fit_fun(data, rand_idx);
    catch
        continue
    end
    residuals = residuals_fun(data, model);
    assert(numel(residuals) == N, 'Residual vector does not match data size')
    inlier_idx = find(residuals < inlier_threshold);
    inlier_ratio = length(inlier_idx)/N;
    if(inlier_ratio > best_inlier_ratio)
        best_model = model;
        best_inlier_ratio = inlier_ratio;
        best_inliers = inlier_idx;
    end
    % Adapt N to inlier ratio
    % See Multiple View Geometry (II ed.) pag.121
    est_outlier_ratio = 1 - inlier_ratio;
    adapt_max_itr = min(idealRansacIterations(0.99, est_outlier_ratio, ...
        sample_size), max_iterations);
    if(inlier_ratio >= inlier_ratio_threshold)
        break
    end
    itr = itr + 1;
end

model = best_model;
inlier_idx = best_inliers;

end

% Prototypes
% model = fit_fun(data, selection)
% residuals = residuals_fun(data, model)
