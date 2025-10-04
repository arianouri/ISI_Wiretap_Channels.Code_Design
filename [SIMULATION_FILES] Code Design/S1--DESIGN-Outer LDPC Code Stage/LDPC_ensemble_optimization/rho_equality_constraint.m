function [eq] = rho_equality_constraint(Rho)
eq = sum(Rho) - 1;