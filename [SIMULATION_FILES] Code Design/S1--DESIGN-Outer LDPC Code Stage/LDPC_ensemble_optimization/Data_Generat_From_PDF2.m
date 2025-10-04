function [ out_random ] = Data_Generat_From_PDF2(pmf_vec,x_vec,num_random)
cdf_vec= cumsum(pmf_vec);
input_random = rand(num_random,1);
out_random = My_lininterpol( x_vec , cdf_vec , input_random );
end

