%% You need to have CVX installed + licensed MOSEK Optimization Tools version 10
% If the reverse counter does not appear in the command window after 5-10 secs
% press Ctrl+C **only a single time** then wait for another 5-10 secs.

clear; clc

load ddist.mat;  %degree distributions

clength     =5e3;

kINP        =1;  %trellis code's input block length
nOUT        =4;  %trellis code's output block length

R_in        =kINP/nOUT;
R_out        =rate_calculation(Lambda,Rho); 

p           =.078;
R_out_p      =R_out/(1-p); %rate of the punctured outer code
R_D         =R_in*R_out_p;  %design rate

R_s         =-0.08+R_D %ultimate secure rate
%                       I(X^N;Z^N)/N = 0.08 bits/cu for this specific ISI-WTC

N_p         =clength/nOUT; %length of the punctured output block

K_s         =ceil(clength*R_s); %length of the secret messages
K_r         =floor(clength*(R_D-R_s)); %length of the dummy messages

N           =round(N_p/(1-p)); %length of the UNpunctured output block

for i=1 %generate multiple H matrices
    
    H = random_PCM_no_girth4(N,Lambda,Rho);
    % [lambda_imp,rho_imp,LAMBDA_imp,GAMMA_imp] = imperical_dd(H);

    save("PCM_"+ num2str(i),'H','-v7.3')
    clear H;
end
