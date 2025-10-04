clear
N_points = 242;
MER_MATRIX = zeros(N_points,3);
save('MER_MATRIX', "MER_MATRIX")
for tstINDex = 1:N_points
    tstINDex
    load(int2str(tstINDex))
    MER_MATRIX(tstINDex,:) = [mean(MAIN_MER) mean(WRTP_MER) mean(WRTP_MER)-mean(MAIN_MER)];
end
    save('MER_MATRIX', "MER_MATRIX")