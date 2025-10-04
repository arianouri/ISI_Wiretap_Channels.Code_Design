function [states, branches_indexed, OUT_CDE]= SupCha_ENCODER(INP_MSG,emMAT,SupCha_CODEBOOK,kINP,nOUT)

    L=length(INP_MSG);
    rs_INP_MSG=(reshape(INP_MSG,kINP,L/kINP)).';
    rs_OUT_CDE=zeros(L/kINP,nOUT);
    

    branches_indexed = zeros(1,L/kINP);
    states = zeros(1,L/kINP);

    currentstate = 1;
    for count = 1:L/kINP
        states(count)=currentstate;
        stage_inp=rs_INP_MSG(count,:);
        
        aux_find_INP_POS=(SupCha_CODEBOOK(:,1)==currentstate)&(SupCha_CODEBOOK(:,[4:4+kINP-1])==stage_inp);
        find_INP_POS=ones(size(aux_find_INP_POS,1),1);
        for andIndex=1:kINP
            find_INP_POS=find_INP_POS & aux_find_INP_POS(:,andIndex);
        end
        INP_POS=find(find_INP_POS==1); 
        rs_OUT_CDE(count,:)=permute(emMAT(SupCha_CODEBOOK(INP_POS,1),SupCha_CODEBOOK(INP_POS,2),SupCha_CODEBOOK(INP_POS,3),:),[4 1 2 3]);
        currentstate=SupCha_CODEBOOK(INP_POS,2);
        branches_indexed(count)=SupCha_CODEBOOK(INP_POS,3);
    end
    OUT_CDE=reshape(rs_OUT_CDE.',1,nOUT*L/kINP);
end