function [true_positive_rate, false_positive_rate, precision, recall] = ...
    detectionStats(problem, inliers)
%detectionStats - Computes some statistics on the results of the detection
%algorithm
%
% Syntax:  [tp, fp, precision, recall] = detectionStats(problem, inliers)
%
% Inputs:
%    problem - A problem structure (see EmptyProblem)
%    inliers - Estimated set of inliers
%
% Outputs:
%    true_positive_rate - True Positive Rate
%    false_positive_rate - False Positive Rate
%    precision - Precision
%    recall - Recall
%
% See: https://en.wikipedia.org/wiki/Precision_and_recall
%
% Author: Pasquale Antonante
% email: antonap@mit.edu
% Date: 2021-01-12

num_outliers = length(problem.outliers);
if(num_outliers == 0)
    true_positive_rate = 1;
    false_positive_rate = 0;
    precision = 1;
    recall = 1;
else
    outliers = setdiff(1:problem.N, inliers);
    true_inliers = setdiff(1:problem.N, union(problem.priors, problem.outliers));
    true_positive = intersect(problem.outliers, outliers);
    false_positive = intersect(true_inliers, outliers);
    true_positive_rate = length(true_positive)/num_outliers;
    false_positive_rate = length(false_positive)/length(true_inliers);
    precision = length(true_positive) / ...
        (length(true_positive) + length(false_positive));
    recall = length(true_positive)/num_outliers;
end
