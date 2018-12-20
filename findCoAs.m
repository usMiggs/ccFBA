%% This script is designed to find CoA containing rxns
%
%   This script will go through a cobra model and identify all rxns
%   containing conenzyme A if it is present on both sides of the rxns. This
%   will avoid flagging the reactions up where coenzyme A gets degraded or
%   synthesiszed. CoA is not contributing to the actual carbon flux and
%   hence should not be taken into account for the carbon constraining
%   algorithm.
%

function [coa] = findCoAs(model,delimiter)

if exist('delimiter','var')
    if isempty(delimiter)
        delimiter = '_';
        warning('default delimiter "_" was chosen!')
    end
else
    delimiter = '_';
    warning('default delimiter "_" was chosen!')
end

% concatenate pattern for regexp
pattern = strcat('coa',delimiter);
pattern2 = '-CoA';
pattern3 = 'coa';

% find all metIDs containing coa
coaCell = regexp(model.mets,pattern);
coaCell2 = regexp(model.metNames,pattern2);
coaCell3 = regexp(model.mets,pattern3);


% convert cell to metID dbl
coaIDs = find(~cellfun(@isempty,coaCell));
coaIDs2 = find(~cellfun(@isempty,coaCell2));
coaIDs3 = find(~cellfun(@isempty,coaCell3));

coaIDs = unique([coaIDs;coaIDs2;coaIDs3]);

% find all rxns that contain a coa metabolite on both sides.
count = 0;
for i = 1:size(model.rxns,1)
    % find substrate and product IDs
    subIDs = find(model.S(:,i) < 0);
    prdIDs = find(model.S(:,i) > 0);
    
    % check if coa met is part of substarte and product list
    subCoA = sum(ismember(subIDs,coaIDs));
    prdCoA = sum(ismember(prdIDs,coaIDs));
    
    % if subCoA and prdCoA are both larger than one it means both sides of
    % the rxn contain the coenzyme and the ID is stored for later
    % processing
    if subCoA > 0 && prdCoA > 0
        count = count + 1;
        coaRxnIDs(count,1) = i;
        
        arr{count,1} = subIDs(ismember(subIDs,coaIDs));
        arr{count,2} = prdIDs(ismember(prdIDs,coaIDs));
        
        subCarb = 0;
        prdCarb = 0;
        
        for j = 1:size(arr{count,1},1)
            subCarb = (abs(full(model.S(arr{count,1}(j,1),i)))*21) + subCarb; % 21 is the carbon per coa
        end
        for j = 1:size(arr{count,2},1)
            prdCarb = (full(model.S(arr{count,2}(j,1),i))*21) + prdCarb; % 21 is the carbon per coa
        end
        
        subCarbon(count,1) = subCarb;
        prdCarbon(count,1) = prdCarb;
        
    end
end

missMatches = find(subCarbon - prdCarbon);



coa.tiles = '*** CoA Rxns ***';
coa.time = date();
coa.demilinter = delimiter;
coa.coaMetIDs = coaIDs;
coa.coaRxnIDs = coaRxnIDs;
coa.subCarbon = subCarbon;
coa.prdCarbon = prdCarbon;
if ~isempty(missMatches)
    coa.missMatch = true;
    coa.numberOfMissMatches = size(missMatches);
    coa.mM_IDs = coaRxnIDs(missMatches);
else
    coa.missMatch = false;
end  












