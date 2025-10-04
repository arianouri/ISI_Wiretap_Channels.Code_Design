clear; clc

load ddist4530_2f3_Rs0427.mat;

clength     =1e5;

R_s         =0.04;

kINP        =1;
nOUT        =4;
N_p         =clength/nOUT;

R_in        =kINP/nOUT;
p           =R_s/(R_in+R_s);

Rout        =rate_calculation(Lambda,Rho);
Rout_p      =Rout/(1-p);

R_D         =R_in*Rout_p;

K_s         =clength*R_s;
K_r         =clength*(R_D-R_s);

N           =N_p+round((K_s*N_p*Rout_p)/(K_s+K_r));

for i=1 %generate multiple H matrices
    
    H = random_PCM_no_girth4(N,Lambda,Rho);
    % [lambda_imp,rho_imp,LAMBDA_imp,GAMMA_imp] = imperical_dd(H);

    save("Rs_0427_2f3_1e5_"+ num2str(i),'H','-v7.3')
    clear H;
end
