function [mtho_pmf]=c_operations_in_SP_or_BP(mthrho,i_pmf,teybel)

i_pmf=i_pmf/sum(i_pmf);
D_c=size(mthrho,2);
mthitho_pmf=i_pmf;
mtho_pmf=zeros(1,size(i_pmf,2));
if mthrho(2)~=0
    mtho_pmf=mthrho(2)*mthitho_pmf;
end
for i=3:D_c,
    mthitho_pmf=basic_c_operations_in_SP_or_BP(i_pmf,mthitho_pmf,teybel);
    if mthrho(i)~=0
        mtho_pmf=mtho_pmf+mthrho(i)*mthitho_pmf;
    end
end