function [ emLLR ] = BCJR_DEC_prevs( SupCha_CODEBOOK , rx_noisless , sigma , Trans_EX , rx_waveform , exLLR , Window_Size)
 State_Cardinality=size(Trans_EX,1);
 numBranches_tot=size(SupCha_CODEBOOK,1);
 nOUT=size(Trans_EX,4);
 kINP=log2(numBranches_tot/State_Cardinality);
 rx_waveform_pac=[reshape(rx_waveform,nOUT,length(rx_waveform)/nOUT)]';
%% Initialization

%extrinsic LLR to probability

CPSideLLR_pac=[reshape(exLLR,kINP,length(exLLR)/kINP)]';
CPSidePr=zeros(size(CPSideLLR_pac,1),kINP,2);
CPSidePr(:,:,1)=(exp(CPSideLLR_pac))./(1+exp(CPSideLLR_pac));
CPSidePr(:,:,2)=1-CPSidePr(:,:,1);



gamma=zeros(numBranches_tot,size(rx_waveform_pac,1));
s_map=zeros(numBranches_tot,3);
%% gamma
ScaleFactor=0;
for state_i=1:State_Cardinality
    for state_j=SpecialiZe(SupCha_CODEBOOK(SupCha_CODEBOOK(:,1)==state_i,2))
        aux_find_bIndex=(SupCha_CODEBOOK(:,[1,2])==[state_i,state_j]);
        find_bIndex=SupCha_CODEBOOK((aux_find_bIndex(:,1) & aux_find_bIndex(:,2))==1,3);
        for branchIndex=find_bIndex'
            ScaleFactor=ScaleFactor+1;
            s_map(ScaleFactor,:)=[state_i,state_j,branchIndex];
%             for X=1:size(rx_noisless,1)
%                 if [permute(Trans_EX(state_i,state_j,branchIndex,:),[4 1 2 3])]'==rx_noisless(X,:)
                      ar=[prod((1./sqrt(2*pi*(sigma^2))).*exp(-(((rx_waveform_pac-[permute(Trans_EX(state_i,state_j,branchIndex,:),[4 1 2 3])]').^2)/(2*(sigma^2)))),2)]';
                      auxPrProd=1;
                      for prodIndex=1:kINP
                          auxPrProd=auxPrProd.*CPSidePr(:,prodIndex,SupCha_CODEBOOK(ScaleFactor,3+prodIndex)+1);
                      end
                      gamma(ScaleFactor,:)=gamma(ScaleFactor,:)+ar.*(auxPrProd');
%                 end
%             end
        end
    end
end

if sum(sum(s_map-SupCha_CODEBOOK(:,[1:3])))~=0
    error("fucked up")
end
% Gtild=log(gamma);

%%
% alpha=zeros(State_Cardinality,size(rx_waveform_pac,1));
% beta=zeros(State_Cardinality,size(rx_waveform_pac,1));
%% alpha (forward recursion)
% alpha(:,1)=1./State_Cardinality;
% alpha(1,1)=1;
% for t=2:size(alpha,2)
%     for state_i=1:State_Cardinality
%         for state_j=SpecialiZe(Trellis_Index_io(Trellis_Index_io(:,1)==state_i,2))
%             alpha(state_j,t)=alpha(state_j,t)+alpha(state_i,t-1).*sum(gamma((s_map(:,1)==state_i)&(s_map(:,2)==state_j),t));
%         end
%     end
% alpha(:,t)=alpha(:,t)/sum(alpha(:,t));
% end
% Atild=log(alpha);
%% beta (backward recursion)
% beta(:,end)=1./State_Cardinality;
% for t=size(beta,2)-1:-1:1
%     for state_i=1:State_Cardinality
%         for state_j=SpecialiZe(Trellis_Index_io(Trellis_Index_io(:,1)==state_i,2))
%             beta(state_j,t)=beta(state_j,t)+beta(state_i,t+1).*sum(gamma((s_map(:,1)==state_i)&(s_map(:,2)==state_j),t+1));
%         end
%     end
% beta(:,t)=beta(:,t)/sum(beta(:,t));
% end
% Btild=log(beta);
%% fast forward backward
post_state_MAR=zeros(numBranches_tot/(2.^kINP),2.^kINP);
past_state_MAR=zeros(size(post_state_MAR));
for iIndex=1:State_Cardinality
    post_state_MAR(iIndex,:)=[SupCha_CODEBOOK(SupCha_CODEBOOK(:,1)==iIndex,2)]'-1;
    past_state_MAR(iIndex,:)=[SupCha_CODEBOOK(SupCha_CODEBOOK(:,2)==iIndex,1)]'-1;
end
past_s_map=zeros(1,numBranches_tot);
nw_sc_index=1;
for auxRindex=1:size(past_state_MAR,1)
    for past_State=SpecialiZe(past_state_MAR(auxRindex,:))
        aux_past_index=[past_State+1,auxRindex]==s_map(:,[1,2]);
        find_aux_past_index=find(aux_past_index(:,1) & aux_past_index(:,2));
        
        past_s_map(nw_sc_index:nw_sc_index+length(find_aux_past_index)-1)=find_aux_past_index-1;
        nw_sc_index=nw_sc_index+length(find_aux_past_index);
    end
end
alpha_fst=zeros(State_Cardinality,Window_Size+1);
beta_fst=zeros(State_Cardinality,Window_Size+1);
alpha_fst(:,1)=1./State_Cardinality;
beta_fst(:,end)=1./State_Cardinality;
[~,~,APP_sigma]=bcjr_fw_bw(alpha_fst,beta_fst,post_state_MAR,past_state_MAR,gamma,Window_Size,int32(past_s_map));
%% APP Pr{i,j|Y}
for index=1:size(APP_sigma,2)
        APP_sigma(:,index)=APP_sigma(:,index)./sum(APP_sigma(:,index));
end
%% Output LLR
APP_pos=zeros(kINP,size(APP_sigma,2));
APP_neg=zeros(size(APP_pos));
for inpIndex=1:kINP
    APP_pos(inpIndex,:)=sum(APP_sigma(SupCha_CODEBOOK(:,3+inpIndex)==0,:));
    APP_neg(inpIndex,:)=sum(APP_sigma(SupCha_CODEBOOK(:,3+inpIndex)==1,:));
end
emLLR=reshape(log(APP_pos./APP_neg),1,length(exLLR));
end