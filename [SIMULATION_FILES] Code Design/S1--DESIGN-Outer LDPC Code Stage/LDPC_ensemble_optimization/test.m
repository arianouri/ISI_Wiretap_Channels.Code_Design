N_points = 1248;
for tstINDex = 1:N_points
    tstINDex
    load(int2str(tstINDex))
    Eq(tstINDex) = ((R_out * R_in)-mean(C_E))/R_s;
end