function [eq] = rate_equality_constraint(lambda)
eq = sum(lambda) - 1;