function problem = EmptyProblem()
%EmptyProblem - Create an empty problem
% Create an empty problem with all the required fields

% Basic elements (these are mandatory)
problem.type= 'none'; % What kind of problem is this?
problem.N = 0; % How many measurements this problem has?
problem.outliers = []; % Among the measurements, what are the indices of the outliers?
problem.priors = [];  % Amont the measurements, what are the indices of the element that will never be outliers?
problem.dof = 0; % Degrees of freedom of the noise process

% Problem specific data
% Add here any problem specific field ...

problem.time = now;
% Unique identifier
problem.uuid = java.util.UUID.randomUUID.toString;
end