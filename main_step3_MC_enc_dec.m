clear all; close all; clc

%% Load Constructed Code

% outer LDPC code stage
load ENC preCalc_p2 inv_T A B permDex H K M N;

% inner trellis code stage
load TRC_42minerr kINP nOUT MAIN WRTP SupCha_CODEBOOK

MAIN = updateSNR(MAIN,-0); % setting the SNR of the main channel = -4

% parameters

R_in    = kINP/nOUT;
R_out   = K/N;

p       =.078;
Rout_p  = R_out/(1-p); %rate of the punctured outer code

R_D     = R_in*Rout_p; %design rate
R_s     =-0.08+R_D %ultimate secure rate
%                   I(X^N;Z^N)/N = 0.08 bits/cu for this specific ISI-WTC
R_r     = R_D-R_s;

Np      = N-round(p*N); %length of the punctured output block
K_s     = floor(nOUT*Np*R_s); %length of the secret messages
K_r     = floor(nOUT*Np*R_r); %length of the dummy messages

K == K_s+K_r % change floor/ceil floor/ceil above if these two were not equal

for sim=1:100 
    iter = 1;
    BER_i = 1;

%% Secret Message

m_s=(rand(1,K_s) < 0.5);

%% Encoding :: outer stage : LDPC code
[cw_outer,cw_outer_punctured] = enc_puncturedLDPC_SC(m_s,preCalc_p2,inv_T,A,B,permDex,p);

punDex = find(isnan(cw_outer_punctured));
unpunDex = find(~isnan(cw_outer_punctured));

cw_outer_punctured_truncated = cw_outer_punctured(unpunDex);

Np == length(cw_outer_punctured_truncated)

%% Encoding :: inner stage : Superchannel
window = 64;
cw_outer_punctured_truncated_bnd = [zeros(1,window+1) cw_outer_punctured_truncated zeros(1,window)];

[~, ~,MAIN.emsSupCha_out]...
    =SupCha_ENCODER(cw_outer_punctured_truncated_bnd,MAIN.SupCha_Trans_EX,SupCha_CODEBOOK,kINP,nOUT);
% [~, ~,WRTP.emsSupCha_out]...
%     =SupCha_ENCODER(cw_outer_punctured_truncated_bnd,WRTP.SupCha_Trans_EX,SupCha_CODEBOOK,kINP,nOUT);

%% Complex Gaussian Noise
MAIN.ChannelOutVec=MAIN.emsSupCha_out...
    +   normrnd(0,MAIN.sigma,1,length(MAIN.emsSupCha_out))...
    +1i*normrnd(0,MAIN.sigma,1,length(MAIN.emsSupCha_out));

% WRTP.ChannelOutVec=WRTP.emsSupCha_out...
%     +   normrnd(0,WRTP.sigma,1,length(WRTP.emsSupCha_out))...
%     +1i*normrnd(0,WRTP.sigma,1,length(WRTP.emsSupCha_out));

%% Turbo Decoding

[oneDex(:,1), oneDex(:,2)]=find(H==1);
E = size(oneDex,1);

MP_MAT = zeros(E,4);                    % 1st column: CN Index
MP_MAT(:,[1 2]) = oneDex;               % 2nd column: VN Index
                                        % 3rd column: V2C mp
                                        % 4rd column: C2V mp
o2i = zeros(1,Np+2*window+1);           % going tb executed by inner
i2o = zeros(1,N);                       % going tb executed by outer
%%

    while  BER_i > 1e-5 && iter < 20
        
        % SC factor graph operations    
            scLLR=BCJR_DEC(SupCha_CODEBOOK,MAIN.sigma,MAIN.SupCha_Trans_EX,MAIN.ChannelOutVec,o2i,window);
            scLLR=reGul(scLLR,33);
        
            i2o(unpunDex) = scLLR(window+2:end-window);
            
            iter
            % disp(['BER inner2outer:',num2str(  sum(sign(i2o)~=1-2*cw_outer)/N  ) ]);
            
        
        for subDex = 1:10
        % VN operations
            [MP_MAT, tx_est] = v_operations(MP_MAT, i2o);
            MP_MAT(:,3)=reGul(MP_MAT(:,3),33);
            
            BER_i = sum(sign(tx_est)~=1-2*cw_outer)/N;
    
            disp(['subIter:', num2str(subDex), ' - BER:',num2str(  BER_i),...
                                               ' - MER:',num2str(  sum(sign(tx_est(punDex))~=1-2*cw_outer(punDex))/K_s  ) ]);
        % CN operation
            MP_MAT = c_operations(MP_MAT, M);
            MP_MAT(:,4)=reGul(MP_MAT(:,4),33);
        end
        
        % VN operations
            [~, vnLLR] = v_operations(MP_MAT, zeros(1,N));
            vnLLR=reGul(vnLLR,33);
        
            o2i(window+2:end-window) = vnLLR(unpunDex);
        
            % i2o(punDex) = vnLLR(punDex);
            
            % disp(['BER outer2inner:',num2str(  sum(sign(vnLLR)~=1-2*cw_outer)/N  ) ]);
            iter = iter+1;
    end
    BER(sim) = sum(sign(tx_est)~=1-2*cw_outer)/N;
    MER(sim) = sum(sign(tx_est(punDex))~=1-2*cw_outer(punDex))/K_s;
end