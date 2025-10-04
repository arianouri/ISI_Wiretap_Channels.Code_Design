function [ Projrction , GapVec ] = myRNG(DataVec,OutLength)
%% this function Generate Random Data Uniformally from Elements of Data Vector
DataVec=unique(DataVec);
randomValues = rand(1,OutLength);
Projrction=zeros(size(randomValues));
for sghIndex=1:OutLength
    Projrction(sghIndex)=...
        DataVec(abs(randomValues(sghIndex)-DataVec)==min(abs(randomValues(sghIndex)-DataVec)));
end
GapVec=abs(Projrction-randomValues);
end

