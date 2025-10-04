function table_gen(maxLLR,number_of_bits,file_name)

deltA=2*maxLLR/power(2,number_of_bits);
a=[deltA:deltA:maxLLR];
tOl=power(2,number_of_bits-1);
teybel=repmat(uint32(0),tOl*(tOl+1)/2,1);
counter=1;
for x=deltA:deltA:maxLLR
    tt=tanh(x/2);
    for y=x:deltA:maxLLR
        tttt=2*atanh(tt*tanh(y/2));
        if tttt >= deltA/2
            if tttt==Inf
                teybel(counter)=min(x,y)/deltA+tOl;
            else
                teybel(counter)=floor(tttt/deltA+0.5)+tOl;  
            end
        else
            teybel(counter)=tOl;
        end
        counter=counter+1;
    end
end
save(file_name,'maxLLR','number_of_bits','teybel')



