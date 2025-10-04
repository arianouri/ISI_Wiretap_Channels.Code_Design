function [ChannelOutVec,NoislessOutVec] = PR_Channel(SourceSymbol,FrequencyResponse,sigma)
% SourceSymbol=-2*(rand(1,500) < 0.5)+1;
SourceSymbol=[1-2*zeros(1,(size(FrequencyResponse,2)-1)),SourceSymbol];
NoislessOutVec=FrequencyResponse(1)*SourceSymbol(size(FrequencyResponse,2):end);
for Index=1:size(FrequencyResponse,2)-1
    NoislessOutVec=NoislessOutVec+FrequencyResponse(Index+1)*SourceSymbol(size(FrequencyResponse,2)-Index:end-Index);
end
% SNR=10^(SNRdB/10);
% sigma=sqrt(sum(FrequencyResponse.^2)./SNR);
% ChannelOutVec= NoislessOutVec + (sigma*randn(1,length(NoislessOutVec)));
% ChannelOutVec=awgn(NoislessOutVec,SNRdB);
ChannelOutVec=NoislessOutVec+normrnd(0,sigma,1,length(NoislessOutVec));
end

