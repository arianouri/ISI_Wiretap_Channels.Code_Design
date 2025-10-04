clear;clc;
load Rule_2;
load TRELLIS_io;
Trans_Pr=TRANS;

numStates=size(Trans_Pr,1);
numBranches=size(Trans_Pr,3);
em_branch = zeros(numStates,length(SpecialiZe(TRELLIS_io(TRELLIS_io(:,1)==1,2)))*numBranches);
kINP=1;
nOUT=4;

%% branch(b_t) pmf, given previous state(s_{t-1})
for iStateIndex=1:numStates
    auxBranchIndex = 1;
    for jStateIndex = [TRELLIS_io(TRELLIS_io(:,1)==iStateIndex,2)]'
        for branchIndex=1:numBranches
            em_branch(iStateIndex,auxBranchIndex) = Trans_Pr(iStateIndex,jStateIndex,branchIndex);
            auxBranchIndex=auxBranchIndex + 1;
        end
    end
end

emBranch_cardinal=size(em_branch,2);
RULE_3_nij_MAT_CUBE=zeros(emBranch_cardinal,numStates,size(RULE_2_suv_kappa_el_MAT,2));
RULE_3_nu_MAT_CUBE=zeros(emBranch_cardinal,numStates,size(RULE_2_suv_kappa_el_MAT,2));
RULE_3_tol_MAT=zeros(numStates,size(RULE_2_suv_kappa_el_MAT,2));

%% OPTIMIZATION
for rl_1_sol_index=1:size(RULE_2_suv_kappa_el_MAT,2) %search among survival answeres of RULE_2
% for rl_1_sol_index=21
    for state_i=1:numStates
    %     for branch_index=1:size(Trans_Pr,3)
            PI_i=(em_branch(state_i,:))';

            cvx_begin
                cvx_quiet true 
            %     cvx_solver sdpt3
            %     cvx_solver SeDuMi
                cvx_solver Mosek
                variable nij(emBranch_cardinal) integer
                minimize( sum(RULE_2_kappa_MAT(:,rl_1_sol_index)...
                    .*sum( kl_div( nij./( 2^kINP*RULE_2_suv_kappa_el_MAT(state_i+1,rl_1_sol_index) ) ,PI_i ) )) )
                sum(nij)==2^kINP*RULE_2_suv_kappa_el_MAT(state_i+1,rl_1_sol_index)
            cvx_end
            RULE_3_tol_MAT(state_i,rl_1_sol_index)=cvx_optval;
    %     end
    RULE_3_nij_MAT_CUBE(:,state_i,rl_1_sol_index)=nij;
    RULE_3_nu_MAT_CUBE(:,state_i,rl_1_sol_index)=nij./( 2^kINP*RULE_2_suv_kappa_el_MAT(state_i+1,rl_1_sol_index) );
    [state_i numStates rl_1_sol_index size(RULE_2_suv_kappa_el_MAT,2)]
    end
end
RULE_3_tol_vec=sum(RULE_3_tol_MAT);

  %  %%  %%  %%  %%
%  %%  %%  %%  %%  %%
  %  %%  %%  %%  %%
    %  %%  %%  %%
      %  %%  %%
        %  %%
          %

optDex = find(RULE_3_tol_vec==min(RULE_3_tol_vec));

RULE_3_AUX_nu_opt=RULE_3_nu_MAT_CUBE(:,:,optDex);
RULE_3_AUX_nij_opt=RULE_3_nij_MAT_CUBE(:,:,optDex);

RULE_3_nu_opt=zeros(size(Trans_Pr));
RULE_3_nij_opt=zeros(size(Trans_Pr));

for iStateIndex=1:numStates
    auxBranchIndex = 1;
    for jStateIndex = [TRELLIS_io(TRELLIS_io(:,1)==iStateIndex,2)]'
        for branchIndex=1:numBranches
            RULE_3_nu_opt(iStateIndex,jStateIndex,branchIndex)=RULE_3_AUX_nu_opt(auxBranchIndex,iStateIndex);
            RULE_3_nij_opt(iStateIndex,jStateIndex,branchIndex)=RULE_3_AUX_nij_opt(auxBranchIndex,iStateIndex);
            auxBranchIndex=auxBranchIndex + 1;
        end
    end
end



%% IMPORTANT: this condition should be satisfied
RULE_3_nij_opt = round(RULE_3_nij_opt);
[sum(RULE_3_nij_opt(1,:,:),'all') sum(RULE_3_nij_opt(:,1,:),'all');...
 sum(RULE_3_nij_opt(2,:,:),'all') sum(RULE_3_nij_opt(:,2,:),'all');...
 sum(RULE_3_nij_opt(3,:,:),'all') sum(RULE_3_nij_opt(:,3,:),'all');...
 sum(RULE_3_nij_opt(4,:,:),'all') sum(RULE_3_nij_opt(:,4,:),'all')]
for tstDex = 1:size(RULE_3_nij_opt,1)
    if   sum(RULE_3_nij_opt(tstDex,:,:),'all')...
       ~=sum(RULE_3_nij_opt(:,tstDex,:),'all')
        error("shit! change 'optDex' and try again")
    end
end



% RULE_3_kappa_opt=RULE_2_suv_kappa_el_MAT(:,optDex);
RULE_3_kappa_opt=round([sum(RULE_3_nij_opt,'all'); sum(RULE_3_nij_opt,[2 3])]/2);

Q=repmat(mu,1,numStates).*Trans_Pr;
Qhat=RULE_3_nij_opt/sum(sum(sum(RULE_3_nij_opt)));
sum(sum(sum(abs(Q-Qhat))))

cd ..
save('RULE_4\Rule_3.mat');