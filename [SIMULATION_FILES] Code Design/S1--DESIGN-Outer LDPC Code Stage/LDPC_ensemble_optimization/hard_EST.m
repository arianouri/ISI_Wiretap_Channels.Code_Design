function [pe] =hard_EST(LLR,TrueSymboles)
HEST=LLR;
HEST(HEST>0)=+1;
HEST(HEST<0)=-1;
pe=1-(sum(HEST==TrueSymboles))/length(TrueSymboles);
end

