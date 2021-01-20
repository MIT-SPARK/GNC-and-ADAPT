function printProblemSummary(problem)
num_outliers = length(problem.outliers);
fprintf('Problem:\n');
fprintf('   Type: `%s`\n', problem.type);
fprintf('   Measurements (N):  %d\n', problem.N);
fprintf('   Num. outliers (k): %d\n', num_outliers);
fprintf('   Num. priors:       %d\n', length(problem.priors));
fprintf('   Num. possible outliers (|V|): %d\n', ...
    problem.N - length(problem.priors));
fprintf('   k/N: %.03f\n', ...
    num_outliers/problem.N);
fprintf('   k/|V|: %.03f\n', ...
    num_outliers/(problem.N - length(problem.priors)));
end