function [ Output ] = Window_BCJR_Decoder(CPSideLLR,RxWaveform,Trans,sigma)
%% Windowed BCJR
%% Variable Initial Size
% Nbit=size(Trans,1)*(2^(-1));
isiLength=log2(size(Trans,1));
CodeLength=size(RxWaveform,2);
WindowSize=(isiLength+1).^2;
CPSidePr_p=zeros(1,CodeLength);
CPSidePr_n=CPSidePr_p;
D=zeros(2^isiLength,2^isiLength,CodeLength);
T=D;
alphaZer=zeros(2^isiLength,CodeLength);
betaZer=alphaZer;
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
betaWin=alphaWin;
Ttr=inv3D(T);
for tIndex=(1+WindowSize):(CodeLength-WindowSize)
    AuxAlphaZer=eye(2^isiLength);
    AuxBetaZer=AuxAlphaZer;
    for sgIndex=1:WindowSize
        AuxAlphaZer=AuxAlphaZer...
                   *(D(:,:,tIndex-sgIndex)*Ttr(:,:,tIndex-sgIndex));
        AuxBetaZer=AuxBetaZer...
                   *inv3D((D(:,:,tIndex+sgIndex)*Ttr(:,:,tIndex+sgIndex)));
    end
    AuxAlphaZer=AuxAlphaZer*alphaWin;
    AuxBetaZer =AuxBetaZer * betaWin;
    alphaZer(:,tIndex)=AuxAlphaZer./sum(AuxAlphaZer);
    betaZer(:,tIndex) =AuxBetaZer./sum(AuxBetaZer);
    
    b(:,tIndex)=(Ttr(:,:,tIndex)*alphaZer(:,tIndex)).*betaZer(:,tIndex);

    b_Num=0;
    b_Den=b_Num;
    
    for bIndex=1:size(TrellisIndexPos+TrellisIndexNeg,2)
       b_Num=b_Num+b(TrellisIndexPos(bIndex),tIndex);
%     end;clear bIndex;
%     for bIndex=1:size(TrellisIndexNeg,2)
       b_Den=b_Den+b(TrellisIndexNeg(bIndex),tIndex);
    end;clear bIndex;
    Output(tIndex)=log(b_Num./b_Den);
end;clear tIndex;
end

