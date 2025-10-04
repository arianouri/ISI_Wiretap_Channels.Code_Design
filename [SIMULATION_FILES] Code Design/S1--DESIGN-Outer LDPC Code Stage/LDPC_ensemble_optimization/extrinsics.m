function [ e2t_pmf_ell ] = extrinsics(Lambda,c2v_pmf_ell)
c2v_pmf_ell=c2v_pmf_ell./sum(c2v_pmf_ell);
i_pmf=c2v_pmf_ell;
tOl=size(c2v_pmf_ell,2);
D_v=size(Lambda,2);
e2t_pmf_ell=zeros(1,tOl);
for j=2:D_v
    i_pmf=basic_v_operations_in_SP_or_BP_or_MS(i_pmf,c2v_pmf_ell);
    if Lambda(j)~=0
        e2t_pmf_ell=e2t_pmf_ell+(Lambda(j)./(j*sum(Lambda./(1:size(Lambda,2))))).*i_pmf;
    end
end
e2t_pmf_ell=e2t_pmf_ell./sum(e2t_pmf_ell);
end

