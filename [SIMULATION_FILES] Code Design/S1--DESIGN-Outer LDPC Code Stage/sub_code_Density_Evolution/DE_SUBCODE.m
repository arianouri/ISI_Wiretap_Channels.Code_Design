function [p_e_in,p_e_out,egzEt_cube_out]=DE_SUBCODE(Rho,AR,Lambda,AL,                                      ...
                                                    sigma,                                                 ...
                                                    SupCha_Trans_EX,SupCha_CODEBOOK,                       ...
                                                    maxLLR,number_of_bits,pE_E,Ell,teybel,                 ...
                                                    TC_inp_length,TC_inp,p,                                ...
                                                    RxWaveform,Window_size)

% AR is an M by 1 vector whose elements identify the algorithms
% 'SP' = sum-product algorithm
% 'BP' = belief propagation algorithm
% 'MS' = min-sum algorithm
%% 'GA' = Gallager's algorithm A
%% 'GB' = Gallager's algorithm B
%% 'MB' = any majority-based algorithm
%% 'AE' = algorithm E

% AL is an N by 1 vector whose elements identify the algorithms
% 'SP' = sum-product algorithm
% 'BP' = belief propagation algorithm
% 'MS' = min-sum algorithm
%% 'GA' = Gallager's algorithm A
%% 'GB' = Gallager's algorithm B
%% 'MB' = any majority-based algorithm
%% 'AE' = algorithm E

bipolar_TC_inp=1-2*TC_inp;

x=linspace(-maxLLR,maxLLR,power(2,number_of_bits)+1);

deltA=2*maxLLR/power(2,number_of_bits);

impulse_zer=zeros(size(x));
impulse_zer(round(maxLLR/deltA)+1)=1;

% impulse_pInf=zeros(size(x));
% impulse_pInf(end)=1;

% impulse_mInf=zeros(size(x));
% impulse_mInf(1)=1;

v2c_pmf_0=impulse_zer;           

egzEt_cube=zeros(Ell,1,length(Lambda));
p_e_1=zeros(1,Ell);
p_e_2=zeros(1,Ell);
p_e = zeros(1,Ell);
p_e_0=sum(v2c_pmf_0(1:round(maxLLR/deltA)))+0.5*v2c_pmf_0(round(maxLLR/deltA)+1);
pe=p_e_0;
ell=1;
while ell<=Ell && pe>=pE_E
    
%% Check-Node to Variable-Node
    if ell==1
        c2v_pmf_ell=c_operations(Rho,AR,v2c_pmf_0,teybel);
    else
        c2v_pmf_ell=c_operations(Rho,AR,v2c_pmf_ell,teybel);
    end
   
%% Variable-Node to Trellis-FG
    v2t_pmf_ell=extrinsics(Lambda,c2v_pmf_ell);
    
%% Trellis-FG to Variable-Node   
    et=Data_Generat_From_PDF2(v2t_pmf_ell,x,TC_inp_length);
    et=bipolar_TC_inp.*et;
    McLLR=BCJR_DEC(SupCha_CODEBOOK,sigma,SupCha_Trans_EX,RxWaveform,et,Window_size);
    McLLR(McLLR>+maxLLR)=+maxLLR;
    McLLR(McLLR<-maxLLR)=-maxLLR;
    McLLR=bipolar_TC_inp.*McLLR;
    
    t2v_pmf_ell=histc(McLLR,x);
    t2v_pmf_ell=t2v_pmf_ell./sum(t2v_pmf_ell);
    
    % uniform puncturing
%     if ell==1
    t2v_pmf_ell=impulse_zer*p+t2v_pmf_ell*(1-p);
%     end

%% Variable-Node to Check-Node
    [v2c_pmf_ell,egzEt_cube(ell,:,:)]=v_operations(Lambda,AL,t2v_pmf_ell,c2v_pmf_ell);
    
%% BER & MER
    Symbol_pmf_ell=basic_v_operations_in_SP_or_BP_or_MS(v2t_pmf_ell,t2v_pmf_ell);
    p_e_1(ell)=sum(Symbol_pmf_ell(1:round(maxLLR/deltA)))+0.5*Symbol_pmf_ell(round(maxLLR/deltA)+1);
    p_e_2(ell)=sum(v2c_pmf_ell(1:round(maxLLR/deltA)))+0.5*v2c_pmf_ell(round(maxLLR/deltA)+1);
    pe=mean([p_e_1(ell),p_e_2(ell)]);
    p_e(ell)=pe;
    [ell p_e_1(ell) p_e_2(ell)]
    ell=ell+1;
end
p_e_in=[p_e_0 p_e(1:sum(p_e~=0)-1)];
p_e_out=p_e(1:sum(p_e~=0));
egzEt_cube_out=egzEt_cube(1:sum(p_e~=0),:,:);
disp(['Number of decoder iterations: ',num2str(ell-1)]);