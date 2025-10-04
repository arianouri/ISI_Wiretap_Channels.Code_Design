clear;clc;
load TRANS;

mu=SState_Pr(TRANS);
S_Cardinal=length(mu);
Kmax=200;
t0L=5e-4;
RULE_2_kappa_el_MAT=zeros(length(mu)+1,Kmax-S_Cardinal+1);
RULE_2_tol_vec=zeros(Kmax-S_Cardinal+1,1);
RULE_2_kappa_MAT=zeros(length(mu),Kmax-S_Cardinal+1);
index=0;

%% OPTIMIZATION
for K=S_Cardinal:Kmax
    index=index+1;
    
    cvx_begin
        cvx_quiet true 
    %     cvx_solver sdpt3
    %     cvx_solver SeDuMi
        cvx_solver Mosek
        variable Ki(S_Cardinal) integer
        minimize( sum( kl_div( Ki/K ,mu ) ) )
        sum(Ki)==K
    cvx_end
    
    RULE_2_tol_vec(index)=cvx_optval;
    RULE_2_kappa_el_MAT(1,index)=K;
    RULE_2_kappa_el_MAT(2:end,index)=Ki;
    RULE_2_kappa_MAT(:,index)=Ki/K;
end
RULE_2_suv_kappa_el_MAT=RULE_2_kappa_el_MAT(:,abs(RULE_2_tol_vec-min(RULE_2_tol_vec))<t0L);
RULE_2_K_opt=RULE_2_kappa_el_MAT(1,abs(RULE_2_tol_vec-min(RULE_2_tol_vec))<t0L);
RULE_2_kappa_opt=RULE_2_kappa_MAT(:,abs(RULE_2_tol_vec-min(RULE_2_tol_vec))<t0L);
cd ..
save('RULE_3\Rule_2.mat');
%           run('RULE_3\Rule_3.m')