function [R] = rate_calculation(Lambda,Rho)
[M,D_c] = size(Rho);
[N,D_v] = size(Lambda);
indices_l = 1./(1:D_v)';
indices_r = 1./(1:D_c)';
R = 1-trace(Rho*repmat(indices_r,1,M))/(Lambda*indices_l);