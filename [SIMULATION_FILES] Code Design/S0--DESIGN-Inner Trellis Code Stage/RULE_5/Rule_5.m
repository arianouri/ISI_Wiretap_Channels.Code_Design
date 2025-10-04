% clear all;clc;
load RULE4_temp;
%% INITIALIZATION
msgLength=kINP*1e5;

%% SOURCE MESSAGE
inp_msg=(rand(1,msgLength) < 0.5);
Window_Size=16;

iterial = 10;
MAIN_MER = zeros(1,iterial);
WRTP_MER = zeros(1,iterial);

for rnDex = 1:iterial
    %% TRELLIS ENCODING
    
    % SupCha_CODEBOOK
    % 1'st col = s_{t-1}
    % 2'nd col = s_t
    % 3'rd col = b_t
    % [4:end]'th col = \{u_0\cdots u_{k-1}\}
    
    NUM_emBranch=2^kINP;
    SupCha_CODEBOOK=[SupCha_trellis(:,[2,3]) ones(size(SupCha_trellis,1),1)];
    for auxIndex=1:size(SupCha_CODEBOOK,1)/2
        emitVec(auxIndex*NUM_emBranch-1:(auxIndex+1)*NUM_emBranch-2)...
            =randperm(2)-1;
    end
    
    SupCha_CODEBOOK=[SupCha_CODEBOOK,emitVec.'];
    

    [~,~,MAIN.emsSupCha_out]...
        =SupCha_ENCODER(inp_msg,MAIN.SupCha_Trans_EX,SupCha_CODEBOOK,kINP,nOUT);
    
    [~,~,WRTP.emsSupCha_out]...
        =SupCha_ENCODER(inp_msg,WRTP.SupCha_Trans_EX,SupCha_CODEBOOK,kINP,nOUT);
    
    %% AWGN CHANNEL
    MAIN.sigma=10.^(-MAIN.SNRdB/20);
    WRTP.sigma=10.^(-WRTP.SNRdB/20);
    
    MAIN.ChannelOutVec=MAIN.emsSupCha_out...
        +   normrnd(0,MAIN.sigma,1,length(MAIN.emsSupCha_out))...
        +1i*normrnd(0,MAIN.sigma,1,length(MAIN.emsSupCha_out));

    WRTP.ChannelOutVec=WRTP.emsSupCha_out...
        +   normrnd(0,WRTP.sigma,1,length(WRTP.emsSupCha_out))...
        +1i*normrnd(0,WRTP.sigma,1,length(WRTP.emsSupCha_out));
    
    %% BCJR DETECTOR
    % isiLength=length(FrequencyResponse)-1;
    % WindowSize=(isiLength+1)^2;
    % FrequencyResponse=[1 -1];
    % isiLength=length(FrequencyResponse)-1;
    % WindowSize=(isiLength+1)^2;
    % [ Trans ] = PR_Trellis( FrequencyResponse );
    % CPSideLLR=zeros(1,length(ChannelOutVec));
    % Output1=fast_Window_BCJR_Detector(CPSideLLR,ChannelOutVec,Trans,sigma,WindowSize);
    
    %% EXTRINSIC INFORMATION
    
    %     inp_msg_LLR(inp_msg==0)=1e2;
    %     inp_msg_LLR(inp_msg==1)=-1e2;
    %     stage_inp_msg_LLR=[reshape(inp_msg_LLR,kINP,length(inp_msg_LLR)/kINP)].';
    %     
    %     stage_CPSideLLR=zeros(size(stage_inp_msg_LLR));
    %     stage_CPSideLLR(:,1)=stage_inp_msg_LLR(:,1);
    %     CPSideLLR=reshape(stage_CPSideLLR.',1,msgLength);
    
        CPSideLLR=zeros(1,msgLength);
    
    %% TRELLIS DECODING
    
    MAIN_Output=BCJR_DEC(SupCha_CODEBOOK,                                   ...
                         MAIN.sigma,                                        ...
                         MAIN.SupCha_Trans_EX,                              ...
                         MAIN.ChannelOutVec,CPSideLLR,Window_Size);
    
    WRTP_Output=BCJR_DEC(SupCha_CODEBOOK,                                   ...
                         WRTP.sigma,                                        ...
                         WRTP.SupCha_Trans_EX,                              ...
                         WRTP.ChannelOutVec,CPSideLLR,Window_Size);
    
    MAIN_stage_Output=[reshape(MAIN_Output,kINP,length(MAIN_Output)/kINP)].';
    
    WRTP_stage_Output=[reshape(WRTP_Output,kINP,length(WRTP_Output)/kINP)].';
    
    
    stage_inp_msg=[reshape(inp_msg,kINP,length(inp_msg)/kINP)].';
    
    % stage_hard_est=(stage_Output<0);
    % MER=sum(stage_hard_est(:,2)~=stage_inp_msg(:,2))/(size(stage_inp_msg,1));
    
    MAIN_hard_est=(MAIN_Output<0);
    WRTP_hard_est=(WRTP_Output<0);
    
    MAIN_MER(rnDex)=sum(MAIN_hard_est~=inp_msg)/(length(inp_msg));
    WRTP_MER(rnDex)=sum(WRTP_hard_est~=inp_msg)/(length(inp_msg));
    
    % save('SUPERCHANNEL.mat','SupCha_CODEBOOK','SupCha_Trans_EX','kINP','nOUT','rx_noisless')
    % save('RULE_5C\Rule_5.mat');

end