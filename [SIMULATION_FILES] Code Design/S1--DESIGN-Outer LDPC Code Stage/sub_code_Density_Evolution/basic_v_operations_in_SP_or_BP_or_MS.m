function [pC]=basic_v_operations_in_SP_or_BP_or_MS(pA,pB)

%pA=pA/sum(pA);
%pB=pB/sum(pB);
tOl=size(pA+pB,2);
L=2*tOl-1;
auxpC=abs(ifft(fft(pA,L).*fft(pB,L),L));
pC=auxpC(round((tOl+1)/2):round((3*tOl-1)/2));
pC(1)=sum(auxpC(1:round((tOl+1)/2)));
pC(tOl)=sum(auxpC(round((3*tOl-1)/2):L));
pC=pC/sum(pC);

