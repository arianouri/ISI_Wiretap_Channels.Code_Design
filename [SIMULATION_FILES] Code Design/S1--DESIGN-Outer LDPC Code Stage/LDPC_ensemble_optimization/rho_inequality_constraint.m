function [ineq, ineq_grad] = rho_inequality_constraint(Rho)
    D_c = length(Rho);
    ineq = zeros(D_c-1,1); 
    ineq_grad = zeros(D_c-1,D_c-1);
    ineq(1:D_c-1) = -Rho(2:end); % positivity of the degree distributions
    ineq_grad((1:D_c-1),:) = -eye(D_c-1);
end