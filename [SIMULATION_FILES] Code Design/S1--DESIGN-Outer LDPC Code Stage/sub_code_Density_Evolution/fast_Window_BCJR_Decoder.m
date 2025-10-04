function [ Output ] = fast_Window_BCJR_Decoder(CPSideLLR,RxWaveform,Trans,sigma,WindowSize)
%% Windowed BCJR
%% Variable Initial Size
% Nbit=size(Trans,1)*(2^(-1));
isiLength=log2(size(Trans,1));
CodeLength=size(RxWaveform,2);
% CPSidePr_p=zeros(1,CodeLength);
% CPSidePr_n=CPSidePr_p;
D=zeros(2^isiLength,2^isiLength,CodeLength);
T=zeros(2^isiLength,2^isiLength,CodeLength);
alphaZer=zeros(2^isiLength,CodeLength);
betaZer=zeros(2^isiLength,CodeLength);
b=zeros(2^isiLength,CodeLength);
Output=zeros(1,CodeLength);
% AuxAlphaZer=eye(2^isiLength);
% AuxBetaZer=AuxAlphaZer;

% Trellis Indexes
sIndex=1;
coIndex=1;
TrellisIndex=zeros(2*2^isiLength,2);
for roIndex=1:2^isiLength
    if coIndex>2^isiLength
        coIndex=1;
    end
    for x=1:2
        TrellisIndex(sIndex,:)=[roIndex,coIndex];
        sIndex=sIndex+1;
        coIndex=coIndex+1;
    end
end

%%
% the state that reach when channel input is -1
    TrellisIndexNeg=SpecialiZe(TrellisIndex(mod(TrellisIndex(:,end),2)==1,end));
% the state that reach when channel input is +1    
    TrellisIndexPos=SpecialiZe(TrellisIndex(mod(TrellisIndex(:,end),2)==0,end));
    CPSidePr_p=(exp(CPSideLLR))./(1+exp(CPSideLLR));
    CPSidePr_n=1-CPSidePr_p;
% for tIndex=1:CodeLength
    for eIndex=1:size(TrellisIndexPos+TrellisIndexNeg,2)
        D(TrellisIndexPos(eIndex),TrellisIndexPos(eIndex),:)=CPSidePr_p;
        D(TrellisIndexNeg(eIndex),TrellisIndexNeg(eIndex),:)=CPSidePr_n;
    end;clear eIndex;   
    for eIndex=1:size(TrellisIndex,1)
        T(TrellisIndex(eIndex,1),TrellisIndex(eIndex,2),:)...
            =exp(-(((RxWaveform-Trans(TrellisIndex(eIndex,1),TrellisIndex(eIndex,2))).^2)./(2*sigma^2)));
    end
% end;clear tIndex;
alphaWin=zeros(2^isiLength,1);
alphaWin(:)=2^(-isiLength);
betaWin=zeros(2^isiLength,1);
betaWin(:)=2^(-isiLength);
sizeTrellisIndexPosNeg=size(TrellisIndexPos+TrellisIndexNeg,2);
[Output] = bcjr_operations(isiLength,...                           %prhs[0]
                           CodeLength,...                          %prhs[1]
                           WindowSize,...                          %prhs[2]
                           D,...                                   %prhs[3]
                           T,...                                   %prhs[4]
                           alphaWin,...                            %prhs[5]
                           betaWin,...                             %prhs[6]
                           alphaZer,...                            %prhs[7]
                           betaZer,...                             %prhs[8]
                           b,...                                   %prhs[9]
                           TrellisIndexPos,...                    %prhs[10]
                           TrellisIndexNeg,...                    %prhs[11]
                           sizeTrellisIndexPosNeg,...             %prhs[12]
                           Output);                               %prhs[13]
Output=log(Output);
Output(1:WindowSize)=0;
Output(CodeLength-WindowSize+1:end)=0;
end

