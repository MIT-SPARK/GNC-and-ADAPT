function n = idealRansacIterations(desired_probability_of_success, outlier_probability, sample_size)
%idealRansacIterations - Computes the optimal number of iterations for
%RANSAC to achieve a desired probability of success
%
% Inputs:
%    desired_probability_of_success - Desired probability of success
%    outlier_probability - Outlier ratio
%    sample_size - The number of samples neede to specify the model
%
% Output:
%    n - the number of iterations
%
% Author: Pasquale Antonante
% email: antonap@mit.edu
% Date: 2021-01-12

assert(desired_probability_of_success > 0, 'The  desired probability  of success must be greater than 0')
n = ceil(log(1-desired_probability_of_success)/log(1-(1-outlier_probability)^sample_size));
assert(isfinite(n), "Infinite number of iterations required")
end

