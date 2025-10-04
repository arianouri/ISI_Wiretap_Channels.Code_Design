function [ntho_pmf,nthegzEt]=v_operations_in_SP_or_BP_or_MS(nthlambda,i_pmf_0,i_pmf)

i_pmf_0=i_pmf_0/sum(i_pmf_0);
i_pmf=i_pmf/sum(i_pmf);
tOl=size(i_pmf_0+i_pmf,2);
D_v=size(nthlambda,2);
nthjtho_pmf=i_pmf_0;
nthegzEt=zeros(1,D_v);
ntho_pmf=zeros(1,tOl);
for j=2:D_v,
    nthjtho_pmf=basic_v_operations_in_SP_or_BP_or_MS(i_pmf,nthjtho_pmf);
    nthegzEt(j)=sum(nthjtho_pmf(1:round((tOl-1)/2)))+0.5*nthjtho_pmf(round((tOl+1)/2));
    if nthlambda(j)~=0
        ntho_pmf=ntho_pmf+nthlambda(j)*nthjtho_pmf;
    end
end
