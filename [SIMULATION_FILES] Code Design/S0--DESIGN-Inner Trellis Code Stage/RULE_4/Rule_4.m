clear;clc;
load Rule_3_60.mat;
load ISIWTC_45_30.mat
%% SUPER-CHANNEL TRELLIS CONSTRUCTION

trns_pr_ht = Qhat./repmat(sum(Qhat,[2 3]),[1 4]);
if min(log2(1./trns_pr_ht),[],'all')<1
    error("adopted transition probabilities fail to suit Kraft-McMillian Ineq.");
end


    SupCha_numEmtBranch=2^kINP;
    SupCha_numStates=RULE_3_kappa_opt(1,1);

    ExCha_numEmtBranch=2^nOUT;
    ExCha_numBranch=size(Trans_Pr,3);
    ExCha_numStates=numStates;

    %%
    % SupCha_trellis is a matrix with columns as follows:

    % 1-Extended channel state at 't-1';
    % 2-Super-channel state at 't-1';
    % 3-Super-channel state at 't';
    % 4-Extended channel state at 't';
    % 5-Branch index '\ell' of extended channel trellis from 't-1' to 't'.
    
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% 

SECAP = 0.107;
for ITERATION=45:50
     
    SupCha_trellis=zeros(SupCha_numStates*SupCha_numEmtBranch,5);
    auxIndex=1;
    for SupCha_state_index=1:SupCha_numStates
        for SupCha_branch_index=1:SupCha_numEmtBranch
            SupCha_trellis(auxIndex,2)=SupCha_state_index;
            auxIndex=auxIndex+1;
        end
    end

    SupCha_state_type=[0;cumsum(RULE_3_kappa_opt(2:end),1)];

    for ExCha_state_index=1:ExCha_numStates
        SupCha_trellis((SupCha_trellis(:,2)>SupCha_state_type(ExCha_state_index))...
                     &(SupCha_trellis(:,2)<=SupCha_state_type(ExCha_state_index+1)),1)=ExCha_state_index;
    end

    %%
    RULE_3_nij_opt=round(RULE_3_nij_opt);
    SupCha_nOrdered_Trellis=SupCha_trellis;
    aux_row_index=1;
    for ex_i_index=1:ExCha_numStates
        posterior=[TRELLIS_io(TRELLIS_io(:,1)==ex_i_index,2)].';
        ex_j_RNDex = randperm(length(posterior));
        for ex_j_index=posterior(ex_j_RNDex)
            for ex_b_index=1:ExCha_numBranch

                if RULE_3_nij_opt(ex_i_index,ex_j_index,ex_b_index)~=0
                    for row_assign=1:RULE_3_nij_opt(ex_i_index,ex_j_index,ex_b_index)
                        SupCha_nOrdered_Trellis(aux_row_index,[4,5])=[ex_j_index ex_b_index];
                        aux_row_index=aux_row_index+1;
                    end
                end

            end
        end
    end
% % TEST#~1
    aux_SupCha_Trellis=SupCha_nOrdered_Trellis(:,[1 4 5]);

    %RULE-4 / REQUIREMENT-3 / RANDOM PERMUTE
    for SupCha_state_index=1:ExCha_numStates
        SupCha_state_index
        aux_SupCha_subTrellis...
            =aux_SupCha_Trellis(...
                SupCha_state_type(SupCha_state_index)...
                    *SupCha_numEmtBranch+1:...
                SupCha_state_type(SupCha_state_index+1)...
                    *SupCha_numEmtBranch,:);
        subTrellis_aux_index=0:SupCha_numEmtBranch:size(aux_SupCha_subTrellis,1);

        % BUF_aux_SupCha_subTrellis=aux_SupCha_subTrellis;
        aux_SupCha_subTrellis = aux_SupCha_subTrellis(randperm(size(aux_SupCha_subTrellis, 1)),:);
        % FLAG=0;
        subTrellis_index=1;
        while subTrellis_index<=length(subTrellis_aux_index)-1
            % if FLAG==1
            %     aux_SupCha_subTrellis=BUF_aux_SupCha_subTrellis;
            % end
            % FLAG=0;
            auxSubTrellis_index = subTrellis_aux_index(subTrellis_index+1);
                overf=0;
                while 1
                    equ_row_index=sim_row(...
                    aux_SupCha_subTrellis(subTrellis_aux_index(subTrellis_index)+1:auxSubTrellis_index,:)...
                   ,aux_SupCha_subTrellis(auxSubTrellis_index,:));

                    if (sum(equ_row_index)==1)
                        break
                    elseif (sum(equ_row_index)==0)&&(~isempty(equ_row_index))
                        error("fucked up");
                    end
                    rnd=round(subTrellis_aux_index(subTrellis_index+1)+1 ...
                        + (size(aux_SupCha_subTrellis,1)...
                        - (subTrellis_aux_index(subTrellis_index+1)+1))*rand);
                    rnd=min([rnd,size(aux_SupCha_subTrellis,1)]);
                    auxBuffer=aux_SupCha_subTrellis(rnd,:);
                    aux_SupCha_subTrellis(rnd,:)=aux_SupCha_subTrellis(auxSubTrellis_index,:);
                    aux_SupCha_subTrellis(auxSubTrellis_index,:)=auxBuffer;
                    overf=overf+1;
                    if (overf>1e4)
                        diffind = sim_row(aux_SupCha_subTrellis(1:auxSubTrellis_index,:)...
                                         ,aux_SupCha_subTrellis(auxSubTrellis_index,:));
                        for lstDex = 1:subTrellis_index-1
                            if diffind(subTrellis_aux_index(lstDex)+1)+diffind(subTrellis_aux_index(lstDex+1))==0
                                lstBuffer=aux_SupCha_subTrellis(subTrellis_aux_index(lstDex)+1,:);
                                aux_SupCha_subTrellis(subTrellis_aux_index(lstDex)+1,:)=aux_SupCha_subTrellis(auxSubTrellis_index,:);
                                aux_SupCha_subTrellis(auxSubTrellis_index,:)=lstBuffer;
                                break
                            end
                        end
                        % FLAG=1;
                        % subTrellis_index=1;
                        % break;
                    end
                end
                % if FLAG==1
                %     break;
                % end
            % end
            % if FLAG==0
                subTrellis_index=subTrellis_index+1;
            % end
        end


    aux_SupCha_Trellis(SupCha_state_type(SupCha_state_index)*SupCha_numEmtBranch+1:SupCha_state_type(SupCha_state_index+1)*SupCha_numEmtBranch,:)=aux_SupCha_subTrellis;
    end
    SupCha_trellis(:,[1 4 5])=aux_SupCha_Trellis;
% % TEST#~2

    % state_vec_Count=SupCha_state_type(1:end-1,1);
    % for rIndex=1:size(SupCha_trellis,1)
    %     if state_vec_Count(SupCha_trellis(rIndex,4))>=SupCha_state_type(SupCha_trellis(rIndex,4)+1)
    %        state_vec_Count(SupCha_trellis(rIndex,4)) =SupCha_state_type(SupCha_trellis(rIndex,4));
    %     end
    %     state_vec_Count(SupCha_trellis(rIndex,4))=state_vec_Count(SupCha_trellis(rIndex,4))+1;
    %     SupCha_trellis(rIndex,3)=R(SupCha_trellis(rIndex,4));    
    % end

for exTrDex = 1:numStates

    st_SC_all = SupCha_trellis(SupCha_trellis(:,1)==exTrDex,2);
    rnd_st_SC_all = st_SC_all(randperm(length(st_SC_all)));

    SupCha_trellis(SupCha_trellis(:,4)==exTrDex,3) = rnd_st_SC_all;

end


% % Test#~3

%% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%% ERROR CHECK    
    errCheck_mat = SupCha_trellis(:,3)==[1:size(SupCha_trellis,1)/2];
    auxfult_VaL=sum(errCheck_mat);
    fult_Index=find(auxfult_VaL~=2^kINP);
    fult_VaL=auxfult_VaL(fult_Index);
    
    downeR=fult_Index(fult_VaL<2^kINP);
    uppeR=fult_Index(fult_VaL>2^kINP);
    
    if sum(sum(SupCha_trellis(:,3)==[1:size(SupCha_trellis,1)/2])~=2^kINP)~=0
        warning('fuckin flood of fuckin shit is comming up')     
        if length(downeR)~=length(uppeR)
            error('flood of shit')
        end
        % 
        % for indeX=1:length(uppeR)
        %     trIndeX=find(SupCha_trellis(:,3)==uppeR(indeX));
        %     SupCha_trellis( (trIndeX(ceil(rand*length(trIndeX)))), 3 )  ...
        %         =downeR(indeX);
        % end
    end
    if sum(sum(SupCha_trellis(:,2)==[1:size(SupCha_trellis,1)/2])~=2^kINP)~=0
        error('flood of shit')
    end
%% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

    K=size(SupCha_trellis,1)/2;
    SupCha_N_numParBranch=1;            
     for SupCha_sIndex=1:K
         for aux_rIndex=1:size(SupCha_trellis(SupCha_trellis(:,2)==SupCha_sIndex,[2,3]),1)
             sub=SupCha_trellis(SupCha_trellis(:,2)==SupCha_sIndex,[2,3]);
            [aux_SupCha_N_numParBranch]=sum(sim_row(sub,sub(aux_rIndex,:)));
            if aux_SupCha_N_numParBranch>SupCha_N_numParBranch
                SupCha_N_numParBranch=aux_SupCha_N_numParBranch;
            end
         end
     end


     MAIN_SupCha_Trans_EX=zeros(K,K,SupCha_N_numParBranch,nOUT);
      MAIN_SupCha_Trans_EX(:)=nan;
       MAIN_Trans_Extended=MAIN.Trans_Extended;

     WRTP_SupCha_Trans_EX=zeros(K,K,SupCha_N_numParBranch,nOUT);
      WRTP_SupCha_Trans_EX(:)=nan;
       WRTP_Trans_Extended=WRTP.Trans_Extended;

     SupCha_Trans_PR=zeros(K,K,SupCha_N_numParBranch);
      SupCha_Trans_PR(:)=nan;

     aux_SupCha=zeros(K,K);

    for SupCha_tr_rIndex=1:size(SupCha_trellis,1)

        aux_SupCha(SupCha_trellis(SupCha_tr_rIndex,2),                      ...
                   SupCha_trellis(SupCha_tr_rIndex,3))                      ...
             =aux_SupCha(SupCha_trellis(SupCha_tr_rIndex,2),                ...
                         SupCha_trellis(SupCha_tr_rIndex,3))+1;

        MAIN_SupCha_Trans_EX( SupCha_trellis(SupCha_tr_rIndex,2) ,          ...
                              SupCha_trellis(SupCha_tr_rIndex,3),           ...
                              aux_SupCha(SupCha_trellis(SupCha_tr_rIndex,2),...
                                SupCha_trellis(SupCha_tr_rIndex,3)),        ...
                              :)                                            ...
                                                                            ...
             =MAIN_Trans_Extended(SupCha_trellis(SupCha_tr_rIndex,1),       ...
                             SupCha_trellis(SupCha_tr_rIndex,4),            ...
                             SupCha_trellis(SupCha_tr_rIndex,5),:);



        WRTP_SupCha_Trans_EX( SupCha_trellis(SupCha_tr_rIndex,2) ,          ...
                              SupCha_trellis(SupCha_tr_rIndex,3),           ...
                              aux_SupCha(SupCha_trellis(SupCha_tr_rIndex,2),...
                                SupCha_trellis(SupCha_tr_rIndex,3)),        ...
                              :)                                            ...
                                                                            ...
             =WRTP_Trans_Extended(SupCha_trellis(SupCha_tr_rIndex,1),       ...
                             SupCha_trellis(SupCha_tr_rIndex,4),            ...
                             SupCha_trellis(SupCha_tr_rIndex,5),:);


        SupCha_Trans_PR( SupCha_trellis(SupCha_tr_rIndex,2) ,               ...
                         SupCha_trellis(SupCha_tr_rIndex,3) ,               ...
                         aux_SupCha(SupCha_trellis(SupCha_tr_rIndex,2),     ...
                         SupCha_trellis(SupCha_tr_rIndex,3)))               ...
                                                                            ...
             =1/size(find(SupCha_trellis(:,2)                               ...
                ==SupCha_trellis(SupCha_tr_rIndex,2)),1);
    end

    MAIN.SupCha_Trans_EX=MAIN_SupCha_Trans_EX;
    WRTP.SupCha_Trans_EX=WRTP_SupCha_Trans_EX;

    SupCha_trellis_io=SupCha_trellis(:,[2,3]);

    %% IUD INFORMATION RATE OF SUPER-CHANNEL
    MONTE_simLength=5e4;
    MONTE_simAvg=20;

    save('MonteCarlo\cube.mat');
    run('MonteCarlo\Main_Of_Mains.m');
%     save(num2str(ITERATION));
 if mean(C_acc)>SECAP
     SECAP = mean(C_acc);
     save('RULE4_temp');
 end
end

cd ..
save('RULE_5\Rule_4.mat');
%           run('RULE_5\Rule_5.m');