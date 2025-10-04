function [rate_obj, rate_obj_grad, rate_obj_hess] = rate_objective(Lambda)
% This function returns the objective function
% Lambda is a row vector, Gradient is returned as a column vector
    D_v = length(Lambda);
    indices_l = 1./(1:D_v)';
    rate_obj = -Lambda*indices_l;
    rate_obj_grad = -indices_l(2:end);
    rate_obj_hess = zeros(D_v-1,D_v-1);
end