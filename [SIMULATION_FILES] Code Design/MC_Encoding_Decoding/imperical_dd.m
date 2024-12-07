function [lambda_imp,rho_imp,LAMBDA_imp,GAMMA_imp] = imperical_dd(H)

H_1=sum(H,1);
LAMBDA_imp = zeros(1,max(H_1));
for i=1:max(H_1)
    LAMBDA_imp(i)=sum(H_1==i)/length(H_1);
end
lambda_imp = ([1:length(LAMBDA_imp)].*LAMBDA_imp)./sum([1:length(LAMBDA_imp)].*LAMBDA_imp);

H_2=sum(H,2);
GAMMA_imp = zeros(1,max(H_2));
for i=1:max(H_2)
    GAMMA_imp(i)=sum(H_2==i)/length(H_2);
end
rho_imp = ([1:length(GAMMA_imp)].*GAMMA_imp)./sum([1:length(GAMMA_imp)].*GAMMA_imp);

end