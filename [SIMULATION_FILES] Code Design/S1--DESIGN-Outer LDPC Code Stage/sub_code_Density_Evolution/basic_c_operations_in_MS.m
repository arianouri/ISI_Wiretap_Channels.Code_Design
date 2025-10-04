function [pC]=basic_c_operations_in_MS(pA,pB)

%pA=pA/sum(pA);
%pB=pB/sum(pB);
tOl=size(pA+pB,2);
PA=cumsum(pA);
AP=fliplr([0 PA(1:tOl-1)]);
PB=cumsum(pB);
BP=fliplr([0 PB(1:tOl-1)]);
PCn=PA.*(1-BP)+(1-AP).*PB;
pCn=PCn-[0 PCn(1:tOl-1)];
pC0=pA(round((tOl+1)/2))+pB(round((tOl+1)/2))-pA(round((tOl+1)/2))*pB(round((tOl+1)/2));
PCp=PA+PB-PA.*PB-AP.*BP;
pCp=PCp-[0 PCp(1:tOl-1)];
pC=[pCn(1:round((tOl-1)/2)) pC0 pCp(round((tOl+3)/2):tOl)];
pC=pC/sum(pC);