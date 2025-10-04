function [ projection ] = Data_Generat_From_PDF1( pdf,RVvalues,OutLength )
%% pdf : the Density function The output Generated Form
%% DataLengh : Number of Symboles
%% OutLengh : Number of Data That you want to be generate correponding to PDF Function
cdf = cumsum(pdf);
cdf=cdf./max(cdf);
% x=linspace(0,DataLength,length(pdf));
[cdf, mask] = unique(cdf);
randomValues=myRNG(cdf,OutLength);
RVvalues=RVvalues(mask);
projection = interp1(cdf, RVvalues, randomValues);
end

