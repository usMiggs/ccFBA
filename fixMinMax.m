function [minmax] = fixMinMax(minmax)
%% fixing the minmax results from runMinMax_GF
%
%   **************************************************
%   *  CREATED by Maximilian Lularevic - 2016-09-16  *
%   **************************************************
%
% this script will scan through a minmax matrix (rxns x 2) and replace all
% non-compliant values with zero (0).
% non-compliant values are values where the min value is larger than the
% max value and vice versa
%
%Example:
%
%   min ---- max
%    0        1
%    1.1      1     <- error
%    -1      -2     <- error
%    0        1
%    0        1



for i = 1:size(minmax,1)
    if minmax(i,1) > minmax(i,2)
        checkBoundsLB(i,1) = 1;
    else
        checkBoundsLB(i,1) = 0;
    end
end

rxnLBissue = find(checkBoundsLB);

for i = 1:size(rxnLBissue,1)
    min = minmax(rxnLBissue(i,1),1);
    max = minmax(rxnLBissue(i,1),2);
    
    minmax(rxnLBissue(i,1),1) = max;
    minmax(rxnLBissue(i,1),2) = min;
end

% for i = 1:size(minmax,1)
%     if minmax(i,1) > minmax(i,2)
%         checkBoundsLB(i,1) = 1;
%     else
%         checkBoundsLB(i,1) = 0;
%     end
% end
% 
% rxnLBissue = find(checkBoundsLB);
% 
% for i = 1:size(rxnLBissue,1)
%     minmax(rxnLBissue(i,1),1) = 0;
% end
% 
% for i = 1:size(minmax,1)
%     if minmax(i,2) < minmax(i,1)
%         checkBoundsUB(i,1) = 1;
%     else
%         checkBoundsUB(i,1) = 0;
%     end
% end    
% 
% rxnUBissue = find(checkBoundsUB);
% 
% for i = 1:size(rxnUBissue,1)
%     minmax(rxnUBissue(i,1),2) = 0;
% end
