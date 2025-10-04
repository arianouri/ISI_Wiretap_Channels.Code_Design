clear all;clc;
Npnt=500;
er_TABLE=zeros(Npnt,2);
for sv_INDEX=243:Npnt
run('Rule_5.m');
[sv_INDEX mean(MAIN_MER) mean(WRTP_MER)]
er_TABLE(sv_INDEX,:)=[mean(MAIN_MER) mean(WRTP_MER)];
save(num2str(sv_INDEX));
end