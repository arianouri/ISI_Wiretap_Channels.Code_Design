function [pC]=basic_c_operations_in_SP_or_BP(pA,pB,teybel)

%pA=pA/sum(pA);
%pB=pB/sum(pB);
tOl=size(pA+pB,2);
pC=check(pA,pB,teybel,tOl);
pC=pC/sum(pC);

