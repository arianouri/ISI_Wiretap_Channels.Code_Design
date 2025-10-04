clear all;clc;

%% 
N=1; %number of the initial points
Input_Size = 3e5;
optIteration=100;

MONTE_simLength=1e5;
MONTE_simAvg=10;

kappa=1;
kappaPrime=3;

wind_length=64;

%%
MAIN.State_Cardinality=2;
WRTP.State_Cardinality=4;
State_Cardinality=max(MAIN.State_Cardinality,WRTP.State_Cardinality);
Number_of_inputs=6; 

%% Signal to Noise Ratio
MAIN.SNRdB=-6.5;
WRTP.SNRdB=-6.0;

%% Channel Frequency Response
% Main Channel
MAIN.FrequencyResponse=[0.792,0.610];
% MAIN.FrequencyResponse=[1 1 -1 -1];
MAIN.FrequencyResponse=MAIN.FrequencyResponse./sqrt(sum(MAIN.FrequencyResponse.^2));
% Eavesdropper's channel
WRTP.FrequencyResponse=[0.445516026180429,0.633021994668546,0.633086585454355];
% WRTP.FrequencyResponse=[1 -1];
WRTP.FrequencyResponse=WRTP.FrequencyResponse./sqrt(sum(WRTP.FrequencyResponse.^2));

%% Joint Source/Wiretap Channel models trellis (EXTENDED)

[ MAIN.Trans_Minimal    , MAIN.Trans_Extended      , MAIN.Trellis_Minimal ,...
  MAIN.Trellis_Extended , MAIN.Trellis_Extended_io , MAIN.branchIndex_Vec ]...
                                                                           ...
     =PR_Trellis_EX(MAIN.FrequencyResponse,Number_of_inputs,MAIN.State_Cardinality);
 

[ WRTP.Trans_Minimal    , WRTP.Trans_Extended      , WRTP.Trellis_Minimal ,...
  WRTP.Trellis_Extended , WRTP.Trellis_Extended_io , WRTP.branchIndex_Vec ]...
                                                                           ...
     =PR_Trellis_EX(WRTP.FrequencyResponse,Number_of_inputs,WRTP.State_Cardinality);
 
 
[ MAIN.Trans_Extended         , WRTP.Trans_Extended,                       ...
  MAIN.Trellis_Extended       , WRTP.Trellis_Extended, TRELLIS,            ...
  MAIN.Trellis_Extended_io    , WRTP.Trellis_Extended_io, TRELLIS_io,      ...
  MAIN.branchIndex_Vec        , WRTP.branchIndex_Vec ]                     ...
                                                                           ...
= ChannelJoinTrellis( MAIN.Trans_Minimal   , WRTP.Trans_Minimal,           ...
                      MAIN.Trellis_Extended , WRTP.Trellis_Extended,       ...
                      MAIN.branchIndex_Vec  , WRTP.branchIndex_Vec );

%% generating equidistributed initial points

% FSMS_S_cardinal=size(MAIN.Trans_Extended,1);

% [p_i_vec]=equidistributed_generator(N,FSMS_S_cardinal);
% [p_i_vec]=repmat(0.5,1,size(MAIN.Trans_Extended,1));


% 
% TRANSVEC=zeros(N,2*FSMS_S_cardinal);
% for iIndex=1:FSMS_S_cardinal
%         TRANSVEC(:,2*iIndex-1)=p_i_vec(:,iIndex);
%         TRANSVEC(:,2*iIndex)=1-p_i_vec(:,iIndex);
% end
% stpCnt=1;
% for MainTR_index=1:N
%     for iIndex=1:2*FSMS_S_cardinal
%             TRANS_0(MainTR_index,TrellisIndex(iIndex,1),TrellisIndex(iIndex,2))=TRANSVEC(MainTR_index,stpCnt);
%             stpCnt=stpCnt+1;
%     end
%     stpCnt=1;
% end


%  TRANS = zeros(State_Cardinality,State_Cardinality,1);
%  for branchIndex=1:size(TRANS,3)
%     TRANS(:,:,branchIndex) = Uniform_Trellis_Generator(size(MAIN.Trans_Extended+WRTP.Trans_Extended,3),TRELLIS_io);
%  end
 
% load TRANS_5th_2
% bufTR=TRANS;
% clear TRANS;
% TRANS(:,:,:)=bufTR(1,:,:,:);
% clear bufTR;


%  TRANS =  [5.77025182116086e-05,0.0129249283102088,0.0803745049546753,0.906642864216904;0.0319402067306963,0.0607722199041158,0.904939802701469,0.00234777066371826;0.00496839049965716,0.401121236624394,0.0603433034958036,0.533567069380145;0.411891608333021,0.584628603114742,0.00240836261914550,0.00107142593309014];
 
TRANS=ones([State_Cardinality,State_Cardinality,size(MAIN.Trans_Extended+WRTP.Trans_Extended,3)]);
TRANS(:,:,:)=TRANS./sum(TRANS,[2,3]);

%% Noiseless Vector

MAIN_Trans_Extended=MAIN.Trans_Extended;
WRTP_Trans_Extended=WRTP.Trans_Extended;

vecIndex=1;
MAIN_rx_noisless_aux=zeros(size(MAIN_Trans_Extended,1)*size(MAIN_Trans_Extended,2)*size(MAIN_Trans_Extended,3),Number_of_inputs);
WRTP_rx_noisless_aux=zeros(size(WRTP_Trans_Extended,1)*size(WRTP_Trans_Extended,2)*size(WRTP_Trans_Extended,3),Number_of_inputs);
for iIndex=1:size(MAIN_Trans_Extended,1)
    for jIndex=1:size(MAIN_Trans_Extended,2)
        for bIndex=1:size(MAIN_Trans_Extended,3)
            MAIN_rx_noisless_aux(vecIndex,:)=[permute(MAIN_Trans_Extended(iIndex,jIndex,bIndex,:),[4 1 2 3])]';
            WRTP_rx_noisless_aux(vecIndex,:)=[permute(WRTP_Trans_Extended(iIndex,jIndex,bIndex,:),[4 1 2 3])]';
            vecIndex=vecIndex+1;
        end
    end
end
clear vecIndex;

MAIN.rx_noisless=unique(MAIN_rx_noisless_aux,'rows');
WRTP.rx_noisless=unique(WRTP_rx_noisless_aux,'rows');

    %% Secrecy Capacity
    
    
%     [TRANSTAR_vec,C1]=Expectation_Maximization_NC_CP( TRANS , ...
%                                           MAIN , WRTP , ...
%                                           Input_Size , kappa , kappaPrime , optIteration , ...
%                                           TRELLIS_io );

    [TRANSTAR_vec,C1]=Expectation_Maximization_MET_i( TRANS , ...
                                          MAIN , WRTP , ...
                                          Input_Size , kappa , kappaPrime , optIteration , ...
                                          TRELLIS_io , wind_length , Number_of_inputs);
%                                       
% 	[TRANSTAR_vec,C1]=Expectation_Maximization_MET_ii( TRANS , ...
%                                           MAIN , WRTP , ...
%                                           Input_Size , kappa , kappaPrime , optIteration , ...
%                                           TRELLIS_io );

C=zeros(optIteration+1,2);  
C(1:end,1)=C1';
mC(1,1)=NaN;

% TRANS=TRANSTAR;
% PERSENT_AGE=0;
save('MonteCarlo\cube.mat');
          run('MonteCarlo\Main_Of_Mains.m');
        
          
% uniq_decodable_bound=min(min(min(log2(1./Trans_Pr_STAR))))/Number_of_inputs;
% save('MonteCarlo\cube.mat');