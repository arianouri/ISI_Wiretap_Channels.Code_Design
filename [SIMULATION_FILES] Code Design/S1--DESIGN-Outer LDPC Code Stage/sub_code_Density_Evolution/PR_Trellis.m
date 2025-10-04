function [ Trans ] = PR_Trellis( FrequencyResponse )
isiLength=size(FrequencyResponse,2)-1;
Nbit=size(FrequencyResponse,2);
tblInput=2*Truth_Table(Nbit)-1;
Trans=zeros(2^isiLength);

sIndex=1;
coIndex=1;
TrellisIndex=zeros(2^Nbit,2);
for roIndex=1:2^isiLength
    if coIndex>2^isiLength
        coIndex=1;
    end
    for x=1:2
%         [roIndex,coIndex]
        Trans(roIndex,coIndex)=sum(inVec(FrequencyResponse).*tblInput(sIndex,:));
        TrellisIndex(sIndex,:)=[roIndex,coIndex];
        sIndex=sIndex+1;
        coIndex=coIndex+1;
    end
end
clear coIndex;
clear roIndex;
end

