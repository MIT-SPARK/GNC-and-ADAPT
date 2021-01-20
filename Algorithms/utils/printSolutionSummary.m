function printSolutionSummary(algname, problem, inliers, info)

[true_positive_rate, false_positive_rate, precision, recall] = ...
    detectionStats(problem, inliers);

outliers = setdiff(1:problem.N, inliers);

fprintf('%s\n', algname);
fprintf('   Iterations: %d\n', info.Iterations);
fprintf('   Correctly estimated outliers: %d / %d\n', ...
    length(intersect(problem.outliers, outliers)), ...
    length(problem.outliers));
if(length(outliers) ~= length(problem.outliers))
    fprintf('   Num. est. outliers: %d\n', length(outliers));
end
fprintf('   True positive rate: %.3f%%\n', ...
    true_positive_rate * 100);
fprintf('   False positive rate: %.3f%%\n', ...
    false_positive_rate * 100);
fprintf('   Precision: %.3f %%\n', precision * 100);
fprintf('   Recall: %.3f %%\n', recall * 100);
if isfield(info, 'time')
    fprintf('   Exec. time: %.2f sec.\n', info.time);
end

end