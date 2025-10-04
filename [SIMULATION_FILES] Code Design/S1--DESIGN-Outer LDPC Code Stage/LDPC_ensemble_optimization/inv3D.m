function [ MDMatT ] = inv3D( MDMat )
MDMatT=zeros(size(MDMat));
for Index=1:size(MDMat,3)
    MDMatT(:,:,Index)=MDMat(:,:,Index)';
end
end

