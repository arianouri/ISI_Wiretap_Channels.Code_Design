function [rate_obj, rate_obj_grad, rate_obj_hess] = rho_objective(Rho)
% This function returns the objective function
% Lambda is a row vector, Gradient is returned as a column vector
    D_c = length(Rho);
    indices_r = 1./(1:D_c)';
    rate_obj = Rho*indices_r;
    rate_obj_grad = indices_r(2:end);
    rate_obj_hess = zeros(D_c-1,D_c-1);
end