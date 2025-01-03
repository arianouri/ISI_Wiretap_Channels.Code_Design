function [cw_outer,cw_outer_punctured] = enc_puncturedLDPC_SC(m_s,preCalc_p2,inv_T,A,B,permDex,p)
   
    K = size(B,2);
    K_s = length(m_s);
    K_r = K-K_s;
    N = length(permDex);


    m_r = (rand(1,K_r) < 0.5);
    s = [m_s, m_r];
    
    p2=mod(preCalc_p2*s',2);
    % P2=;
    p1=mod(-inv_T*(mod(A*p2,2)+mod(B*s',2)),2);
    
    cw_outer_sys=[p1; p2; s'];
    
    m_s=m_s+0;
    m_s(randperm(K_s,round(p*N)))=NaN;
    cw_outer_punctured_sys=[p1; p2; [m_s m_r]'];
    
    [~,truDex]=sort(permDex);

    cw_outer(1:N)=cw_outer_sys(truDex);
    cw_outer_punctured(1:N)=cw_outer_punctured_sys(truDex);
end