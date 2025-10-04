load 1.mat

mean_MAIN_p_e_in = zeros(1,300);
mean_MAIN_p_e_in(1:length(MAIN_p_e_in)) = MAIN_p_e_in;

mean_FICT_p_e_in = zeros(1,300);
mean_FICT_p_e_in(1:length(FICT_p_e_in)) = FICT_p_e_in;

for inDex = 2:5
    load(num2str(inDex))

    mean_MAIN_p_e_in = mean([mean_MAIN_p_e_in; MAIN_p_e_in zeros(1,300-length(MAIN_p_e_in))],1);
    mean_FICT_p_e_in = mean([mean_FICT_p_e_in; FICT_p_e_in zeros(1,300-length(FICT_p_e_in))],1);
end
mean_FICT_p_e_in = [mean_FICT_p_e_in mean_MAIN_p_e_in(end)];

plot(mean_MAIN_p_e_in(mean_MAIN_p_e_in>eps)); hold on
plot(mean_FICT_p_e_in(mean_FICT_p_e_in>eps));
plot(0.5*ones(size(mean_MAIN_p_e_in)));
