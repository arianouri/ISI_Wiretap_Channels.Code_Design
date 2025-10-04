clc;
clear all;
load('42minerr.mat')
save('SUPERCHANNEL.mat','kINP','nOUT','MAIN','WRTP','SupCha_CODEBOOK')
save('RateLimits.mat','C_acc','C_B','C_E')
clear

load('SUPERCHANNEL.mat');
load('RateLimits.mat');
TC_inp_length=2.5e4;

Window_size=64;


%% DE parameters

% maximal log likelihood ratio (LLR)
maxLLR=30;
% number of bits for quantizing the interval [-maxLLR:maxLLR]
number_of_bits=13;
% name of the file in which tanh table 'teybel' is saved
file_name='lookup_table';
% tanh table 'teybel' is generated
table_gen(maxLLR,number_of_bits,file_name);
% tanh table 'teybel' is loaded into memory
load('lookup_table');

% error floor
stop_pe=10^-5 ;
% maximal number of iterations
iter=300;

%% tolerance in the error of the optimization function
tol = 5e-3;  
tol_feas = 5e-3;

%% Rates

% C_s               =mean(C_acc);




% p                 =
% R_s               =R_in*p/(1-p);

R_SC_s            =0.11;
R_in              =kINP/nOUT;

R_s               =R_SC_s;
p                 =R_s/(R_in+R_s);

R_SC_BOB          =mean(C_B);
% R_SC_BOB           =0.101;
% R_SC_EVE           =R_SC_BOB-R_SC_s;

R_out_target      =R_SC_BOB/R_in;

%% Degree Distributions

load ddist4530_453_c2_Rs07;

% load ddistPEG;
% Lambda = lambda_imp;
% Rho = rho_imp;

% algorithms are identified
AR=['SP'];   % or 'BP' for sum-product or belief propagation algorithms];     

%% PHASE-I
% 1 by D_v variable-nodes' distribution, D_v is the maximal variabe degree
% RR=R_out;
% sumRho=sum(Rho./(1:length(Rho)));
% D_v=40; %maximum degree of variable node
% Lrow=1./(1:D_v);
% cvx_begin
%     cvx_quiet true 
% %     cvx_solver sdpt3 
% %     cvx_solver SeDuMi
%     cvx_solver Mosek
%     variable Lambda(D_v)
%     minimize (norm((1-RR)*Lrow*Lambda-sumRho))
%     subject to
%         Lambda >= 0
%         sum(Lambda) == 1
%         Lambda(1) == 0
% %         Lambda(D_v) == 0.2
% cvx_end
% Lambda=Lambda';
% if rate_calculation(Lambda,Rho)-RR>1e-6 
%     error('PHASE-I error'); end

% algorithms are identified
AL=['SP']; % or 'BP' or 'MS' for sum-product or belief propagation or min-sum algorithms

%%

R_out = rate_calculation(Lambda,Rho);
R_out = R_out/(1-p);

disp(['Design Rate: ',num2str(  (R_out * R_in)  ) ]);
disp(['Equivocation Rate: ',num2str(  ((R_out * R_in)-mean(C_E))/R_s  )]);

for simDex = 2:100

%% Monte-Carlo
%% Trellis Encoding
% Trellis input / LDPC output
TC_inp_msg=(rand(1,TC_inp_length) < 0.5);
% Puncturing
% TC_inp_msg_PUNCTURED=[(rand(1,k_s) < 0.5), TC_inp_msg(k_s+1:end)];
% encoder
[~, ~,MAIN.emsSupCha_out]...
    =SupCha_ENCODER(TC_inp_msg,MAIN.SupCha_Trans_EX,SupCha_CODEBOOK,kINP,nOUT);
[~, ~,WRTP.emsSupCha_out]...
    =SupCha_ENCODER(TC_inp_msg,WRTP.SupCha_Trans_EX,SupCha_CODEBOOK,kINP,nOUT);

%% Channel Simulation

MAIN.ChannelOutVec=MAIN.emsSupCha_out...
    +   normrnd(0,MAIN.sigma,1,length(MAIN.emsSupCha_out))...
    +1i*normrnd(0,MAIN.sigma,1,length(MAIN.emsSupCha_out));

WRTP.ChannelOutVec=WRTP.emsSupCha_out...
    +   normrnd(0,WRTP.sigma,1,length(WRTP.emsSupCha_out))...
    +1i*normrnd(0,WRTP.sigma,1,length(WRTP.emsSupCha_out));

%%
%     Rho = LDPC_opt_RHO('rho_objective',                                 ...
%                    'rho_inequality_constraint',                         ...
%                    'rho_equality_constraint',                           ...
%                     Rho,Lambda,100);

% Density Evolution %% % % % % % % % % % % % % % % % % % % % % % % % % % % 

[MAIN_p_e_in,~,~]         =DE_SUBCODE(Rho,AR,Lambda,AL,                                ...
                                      MAIN.sigma,                                       ...
                                      MAIN.SupCha_Trans_EX,SupCha_CODEBOOK,             ...                                 ...
                                      maxLLR,number_of_bits,stop_pe,iter,teybel,        ...
                                      TC_inp_length,TC_inp_msg,p,                       ...
                                      MAIN.ChannelOutVec,Window_size)
                                      
% [WRTP.p_e_in,~,~]         =DE_SUBCODE(Rho,AR,Lambda,AL,                                 ...
%                                       WRTP.sigma,                                       ...
%                                       WRTP.SupCha_Trans_EX,SupCha_CODEBOOK,             ...
%                                       maxLLR,number_of_bits,stop_pe,iter,teybel,        ...
%                                       TC_inp_length,TC_inp_msg,p,                       ...
%                                       WRTP.ChannelOutVec,Window_size);
                                  
[FICT_p_e_in,~,~]    =DE_SUBCODE_FICT(Rho,AR,Lambda,AL,                                 ...
                                      WRTP.sigma,                                       ...
                                      WRTP.SupCha_Trans_EX,SupCha_CODEBOOK,             ...
                                      maxLLR,number_of_bits,stop_pe,iter,teybel,        ...
                                      TC_inp_length,TC_inp_msg,p,                       ...
                                      WRTP.ChannelOutVec,Window_size);

save(num2str(simDex));
end