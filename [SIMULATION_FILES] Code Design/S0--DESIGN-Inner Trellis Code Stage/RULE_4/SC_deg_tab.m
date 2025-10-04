clear all; clc;
load('0044962_itr100.mat');
sp_deg_tab=zeros(size(Q,1)*size(Q,2)*size(Q,3),3);
fl_counter=0;
for i=1:size(Q,1)
    for k=1:size(Q,2)
        for j=1:size(Q,3)
            fl_counter=fl_counter+1;
            sp_deg_tab(fl_counter,:)=[Q(i,j,k) RULE_3_nij_opt(i,j,k) Qhat(i,j,k)];
        end
    end
end
sp_deg_tab_rounded=round(sp_deg_tab, 3);