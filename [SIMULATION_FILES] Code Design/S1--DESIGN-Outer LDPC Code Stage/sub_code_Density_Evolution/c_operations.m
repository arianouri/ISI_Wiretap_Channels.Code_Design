function [o_pmf]=c_operations(Rho,AR,i_pmf,teybel)

% Rho =
%       [0      rho2_1  rho3_1  ...     rhoD_c_1
%        0      rho2_2  rho3_2  ...     rhoD_c_2
%        .      .       .       ...     .
%        .      .       .       ...     .
%        .      .       .       ...     .
%        0      rho2_M  rho3_M  ...     rhoD_c_M]

% AR is an M by 1 vector whose elements identify the algorithms
% 'SP' = sum-product algorithm
% 'BP' = belief propagation algorithm
% 'MS' = min-sum algorithm
%% 'GA' = Gallager's algorithm A
%% 'GB' = Gallager's algorithm B
%% 'MB' = any majority-based algorithm
%% 'AE' = algorithm E

i_pmf=i_pmf/sum(i_pmf);
M=size(Rho,1);
o_pmf=zeros(1,size(i_pmf,2));
for m=1:M,
    if AR(m,:)=='SP' | AR(m,:)=='BP'
        mtho_pmf=c_operations_in_SP_or_BP(Rho(m,:),i_pmf,teybel);
    elseif AR(m,:)=='MS'
        mtho_pmf=c_operations_in_MS(Rho(m,:),i_pmf);
%%   elseif AR(m,:)=='GA' | AR(m,:)=='GB' | AR(m,:)=='MB'
%%       MAPPER
%%       mtho_pmf=c_operations_in_GA_or_GB_or_MB(Rho(m,:),i_pmf);
%%       DEMAPPER
%%   elseif AR(m,:)=='AE'
%%       MAPPER
%%       mtho_pmf=c_operations_in_AE(Rho(m,:),i_pmf);
%%       DEMAPPER
    else
        warning([AR(m,:) '! what on earth is ' AR(m,:) '?']);
        MM;
    end
    o_pmf=o_pmf+mtho_pmf;
end
o_pmf=o_pmf/sum(o_pmf);

