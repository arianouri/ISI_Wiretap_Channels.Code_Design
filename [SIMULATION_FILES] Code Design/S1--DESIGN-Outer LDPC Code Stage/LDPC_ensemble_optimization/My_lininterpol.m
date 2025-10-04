function out = My_lininterpol( x_vec , cdf_vec , in_vec )
out=zeros(1,length(in_vec));
    for i=1:length(in_vec)
        if in_vec(i) <= cdf_vec(1)
            out(i) = x_vec(1) ;

        elseif in_vec(i) >= cdf_vec(end)
            out(i) = x_vec(end) ;

        else
%             Indtrash=[];Ind=[];
%             Indtrash = find(  cdf_vec < in_vec(i) ) ;
%             Ind = max(Indtrash);
            Ind = find(  cdf_vec < in_vec(i) ,1, 'last' ) ;
           
            m = ( x_vec(Ind+1) - x_vec(Ind) ) / ( cdf_vec(Ind+1) - cdf_vec(Ind) ) ;
            out(i) = x_vec(Ind) + m*( in_vec(i) - cdf_vec(Ind) ) ;
        end
    end
end