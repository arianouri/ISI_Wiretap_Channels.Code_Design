function [ DEC ] = LutToDec( TABLE )
TABLE=(TABLE+1)./2;
for Index=1:size(TABLE,1)
    StringAux=num2str(TABLE(Index,:));
    StringAux(StringAux==' ')='';
    String(Index,:)=StringAux;
end
DEC=bin2dec(String);
end

