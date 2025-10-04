function [o_pmf,egzEt]=v_operations(Lambda,AL,i_pmf_0,i_pmf)

% Lambda =
%          [0       lambda2_1   lambda3_1   ...     lambdaD_v_1
%           0       lambda2_2   lambda3_2   ...     lambdaD_v_2
%           .       .           .           ...     .
%           .       .           .           ...     .
%           .       .           .           ...     .
%           0       lambda2_N   lambda3_N   ...     lambdaD_v_N]

% AL is an N by 1 vector whose elements identify the algorithms
% 'SP' = sum-product algorithm
% 'BP' = belief propagation algorithm
% 'MS' = min-sum algorithm
%% 'GA' = Gallager's algorithm A
%% 'GB' = Gallager's algorithm B
%% 'MB' = any majority-based algorithm
%% 'AE' = algorithm E

i_pmf_0=i_pmf_0/sum(i_pmf_0);
i_pmf=i_pmf/sum(i_pmf);
tOl=size(i_pmf_0+i_pmf,2);
[N,D_v]=size(Lambda);
egzEt=zeros(N,D_v);
o_pmf=zeros(1,tOl);
for n=1:N,
    if AL(n,:)=='SP' | AL(n,:)=='BP' | AL(n,:)=='MS'
        [ntho_pmf,nthegzEt]=v_operations_in_SP_or_BP_or_MS(Lambda(n,:),i_pmf_0,i_pmf);
%%   elseif AL(n,:)=='GA' | AL(n,:)=='GB' | AL(n,:)=='MB'
%%       MAPPER
%%       [ntho_pmf,nthegzEt]=v_operations_in_GA_or_GB_or_MB(Lambda(n,:),i_pmf_0,i_pmf);
%%       DEMAPPER
%%   elseif AL(n,:)=='AE'
%%       MAPPER
%%       [ntho_pmf,nthegzEt]=v_operations_in_AE(Lambda(n,:),i_pmf_0,i_pmf);
%%       DEMAPPER
    else
        warning([AL(n,:) '! what on earth is ' AL(n,:) '?']);
        NN;
    end
    o_pmf=o_pmf+ntho_pmf;
    egzEt(n,:)=nthegzEt;
end
o_pmf=o_pmf/sum(o_pmf);

