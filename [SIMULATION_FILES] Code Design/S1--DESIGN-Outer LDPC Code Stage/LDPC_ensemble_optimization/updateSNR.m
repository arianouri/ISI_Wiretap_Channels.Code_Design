function [STR] = updateSNR(STR,snr)

    STR.SNRdB = snr;
    STR.sigma = 10^(-(snr)/20);
    STR.VarEyAns = (STR.sigma)^2;

end