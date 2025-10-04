%% TEST#~1

% for cIndex = 1:4
%     for rIndex = 1:4
%         for bIndex = 1:4
% % [a b]
%             a = length(find(SupCha_nOrdered_Trellis(:,1)==cIndex & SupCha_nOrdered_Trellis(:,4)==rIndex & SupCha_nOrdered_Trellis(:,5)==bIndex));
%             b = RULE_3_nij_opt(cIndex,rIndex,bIndex);
%             if a~=b
%                 [a b]
%             end
% 
%         end
%     end
% end

%% Test#~2

% for cIndex = 1:4
%     for rIndex = 1:4
%         for bIndex = 1:4
% [a b]
%             a = length(find(SupCha_trellis(:,1)==cIndex & SupCha_trellis(:,4)==rIndex & SupCha_trellis(:,5)==bIndex));
%             b = RULE_3_nij_opt(cIndex,rIndex,bIndex);
%             if a~=b
%                 [a b]
%             end
% 
%         end
%     end
% end

%% Test#~3
% 
% for Dex=1:200
%     fnd = find(SupCha_trellis(:,3)==Dex);
%     if length(fnd)~=2
%         error("damn")
%     end
% end
