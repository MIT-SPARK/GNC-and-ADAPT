function problem = linearRegressionProblem(N, outliers_percentage, varargin)
% Generates a random linear estimation problem (e.g. like ||Ax-y||)
% 
% Inputs:
% - N: number of measurements
%
% Options:
% - 'StateDimension': size of x (default: 3)
% - 'MeasurementDimension': size of y (default: 3)
% - 'StateMeanAndStd': x mean and std. (default: [1, 1])
% - 'MeasurementNoiseStd': y std. (default: 0.1)
% - 'OutliersSpace': range for an outlier
%
% Outputs:
% - The estimation problem including
%   - type: 'linear-estimation'
%   - N: number of measurements
%   - StateDimension: dimension of the variable to be estimated
%   - MeasurementDimension: dimension of each measurement
%   - x_gt: the ground truth state
%   - A: cell array containing each A, (i.e. y_i = A_i*x)
%   - y: cell array containing each measurements plus noise
%   - outliers: indices of outlier measurements
%   - priors: empty vector
%   - uuid: unique identifier
%
%
% Edited: Pasquale Antonante
% Date: 2021-01-12
% Institution: Massachusetts Institute of Technology

%% Parameters
params = inputParser;
params.CaseSensitive = false;
params.addParameter('StateDimension', 3, ...
    @(x) isnumeric(x) && isscalar(x) && x>=1);
params.addParameter('MeasurementDimension', 3, ...
    @(x) isnumeric(x) && isscalar(x) && x>=1);
% Position (x-y) parameters
params.addParameter('StateMeanAndStd', [1 1], ...
    @(x) isnumeric(x) && isvector(x) && numel(x) == 2 && all(x>=0));
params.addParameter('MeasurementNoiseStd', 0.1, ...
    @(x) isnumeric(x) && isscalar(x) && x>0);
params.addParameter('OutliersSpace', [-100 100], ...
    @(x) isnumeric(x) && isvector(x) && numel(x) == 2 && x(1) < x(2));
% Parse 'em all!
params.parse(varargin{:});

% aliases
x_dim = params.Results.StateDimension;
y_dim = params.Results.MeasurementDimension;
x_mean_std = params.Results.StateMeanAndStd;
noise_std = params.Results.MeasurementNoiseStd;
outliers_space = params.Results.OutliersSpace;

if nargin < 2
    outliers_percentage = 0;
end

%% Generation
problem = EmptyProblem();

problem.N = N; % num samples
problem.StateDimension = x_dim;
problem.MeasurementDimension = y_dim;
problem.MeasurementNoiseStd = noise_std;
problem.dof = y_dim;
problem.type = 'linear-estimation';

num_outliers = floor(outliers_percentage*N);
% x vector, randomly generated
problem.x_gt = normrnd(x_mean_std(1), x_mean_std(2), x_dim, 1);

% A Matrix, randomly generated for each y
problem.A = {};
for i=1:problem.N
    problem.A{i} = random_matrix(y_dim, x_dim, min([x_dim y_dim]));
end

% y Measurements
problem.y = {};
for i=1:problem.N
    problem.y{i} = problem.A{i}*problem.x_gt + normrnd(0, noise_std, y_dim, 1);
end

% Add outliers
convi = 1/noise_std*eye(y_dim);
problem.outliers = sort(randperm(problem.N, num_outliers));
for i = problem.outliers  
    y_gt = problem.y{i};
    d = 0;
    while(d < chi2inv(0.999, y_dim))
        y_meas = rand_in(y_dim, 1, outliers_space);
        d = (y_gt - y_meas)'*convi*(y_gt - y_meas);
    end
    problem.y{i} = y_meas;
end

% Priors
problem.priors = [];

end

%% Utilities
function r = rand_in(n, m, interval)
    a = interval(1);
    b = interval(2);
    r = a + (b-a).*rand(n,m);
end

function r = randn_in(n, m, interval)
    a = interval(1);
    b = interval(2);
    r = a + (b-a).*randn(n,m);
end

function M = random_matrix(n, m, desired_rank)
    max_rank = min([n m]);
    if nargin<3
        desired_rank = max_rank;
    end
    assert(desired_rank <= max_rank, ...
        'Cannot generate a matrix a random matrix with the desired rank')

    M = 3*rand(n, m);
    while rank(M) < max_rank && rank(M) > 0 % need to start from a full-rank matrix
        M = randn_in(n, m, [-1, 1]);
    end

    % decrease rank if needed
    rank_M = rank(M); % this is just and estimation of the rank
    if rank_M ~= desired_rank
        col_perm = randperm(m, rank_M - desired_rank  + 1);
        for i=2:(rank_M - desired_rank  + 1)
            M(:,col_perm(1)) = M(:,col_perm(i));
        end
    end
    
end