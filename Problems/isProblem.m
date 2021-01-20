function flag = isProblem(problem)
% Check the input problem has all the mandatory fields
mandatory_fields = {'type', 'N', 'outliers', 'priors', 'dof'};
flag = all(isfield(problem, mandatory_fields));
end

