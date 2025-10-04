clear all;close all;clc;
load cube.mat
collector_samples=0;
Input_Size = MONTE_simLength ;  
MeanIter = MONTE_simAvg;
% Trans_Pr=zeros(size(Trans_Extended,1),size(Trans_Extended,2),size(Trans_Extended,3));

    %% MAin Channel Parameters

    MAIN.VarEyAns=10^(-(MAIN.SNRdB)/10);
    WRTP.VarEyAns=10^(-(WRTP.SNRdB)/10);

    %%
    format long
    C_acc=zeros(1,MeanIter);
    C_B=zeros(1,MeanIter);
    C_E=zeros(1,MeanIter);
        for mIndex=1:MeanIter
                
            C_B(mIndex)=binary_rate( Input_Size , MAIN.VarEyAns , MAIN.SupCha_Trans_EX , SupCha_Trans_PR , SupCha_trellis_io );
            C_E(mIndex)=binary_rate( Input_Size , WRTP.VarEyAns , WRTP.SupCha_Trans_EX , SupCha_Trans_PR , SupCha_trellis_io );
            
            C_acc(mIndex)=C_B(mIndex)-C_E(mIndex);
            
            disp(['Capacity Sim','(',num2str(mIndex),') :',num2str(C_acc(mIndex))]);
            symbline(length(['Capacity Sim','(',mIndex,') :',num2str(C_acc(mIndex))]),'-');
        
        end
   
   disp(['Eb/N0 BOB: ',num2str(MAIN.SNRdB),'Eb/N0 EVE: ',num2str(WRTP.SNRdB),' | Capacity Avg: ',num2str(mean(C_acc))]);
   symbline(length(['Eb/N0 BOB: ',num2str(MAIN.SNRdB),'Eb/N0 EVE: ',num2str(WRTP.SNRdB),' | Capacity Avg: ',num2str(mean(C_acc))]),'=');
