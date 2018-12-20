%% This script is designed to find large molecule containing rxns
%
%   This script will go through a cobra model and identify all rxns
%   containing large molecules meaning cytochromes, hemes, and quinones
%   which do not contriobute to the carbon flux of a reaction but rather as
%   an electron carrier.
%   Only rxns where these molecules are present on both sides of the rxns, 
%   i.e. their oxidized and reduced will be identified. This will avoid 
%   flagging the reactions up where they get degraded or synthesiszed. 
%

function [lrgMlcs] = findLrgMlcs(model,carbonCount)

if exist('carbonCount','var')
    if isempty(carbonCount)
        carbonCount = findMetCarbon(model);
    elseif size(carbonCount.carbonAtoms,1) < size(model.mets,1) || size(carbonCount.carbonAtoms,1) > size(model.mets,1)
        carbonCount = findMetCarbon(model);
    end
else
    carbonCount = findMetCarbon(model);
end

% find all metIDs containing large molecules
quiCell = regexp(model.metNames,'quino');
cytCell = regexp(model.metNames,'cytoch');
hemCell = regexp(model.metNames,'heme');

% convert cell to metID dbl
quiIDs = find(~cellfun(@isempty,quiCell));
cytIDs = find(~cellfun(@isempty,cytCell));
hemIDs = find(~cellfun(@isempty,hemCell));

lrgIDs = [quiIDs;cytIDs;hemIDs];

rxnIDs = [];
for i = 1:size(lrgIDs,1)
    tmp = [];
    tmp = find(model.S(lrgIDs(i,1),:));
    rxnIDs = [rxnIDs;tmp'];
end

[C,IA,IC] = unique(rxnIDs);
valCount = hist(rxnIDs(:,1),C)';
trgtRxns = C(valCount>1);

% find logical vector for transportRxns
isTrans  = findTrans(model);

%retrieve actual reactions containing electron transferring mlcls
lrgMlcRxns= trgtRxns(~ismember(trgtRxns,find(isTrans)));

for i = 1:size(lrgMlcRxns,1)
    % find substrate and product IDs
    subIDs = find(model.S(:,lrgMlcRxns(i,1)) < 0);
    prdIDs = find(model.S(:,lrgMlcRxns(i,1)) > 0);
    
    % check if coa met is part of substarte and product list
    
    subLM = subIDs(ismember(subIDs,lrgIDs));
    prdLM = prdIDs(ismember(prdIDs,lrgIDs));
    
    if size(subLM,1) > 1 && size(prdLM,1) > 1
        for j = 1:size(subLM,1)
            subCarb(j,1) = carbonCount.carbonAtoms(subLM(j,1),1)*...
                abs(full(model.S(subLM(j,1),lrgMlcRxns(i,1))));

            prdCarb(j,1) = carbonCount.carbonAtoms(prdLM(j,1),1)*...
                abs(full(model.S(prdLM(j,1),lrgMlcRxns(i,1))));
        end

        subCarb = sum(subCarb);
        prdCarb = sum(prdCarb);

    elseif size(subLM,1) > size(prdLM,1)
        if size(prdLM,1) > 1
            for j = 1:size(subLM,1)
                subCarb(j,1) = carbonCount.carbonAtoms(subLM(j,1),1)*...
                    abs(full(model.S(subLM(j,1),lrgMlcRxns(i,1))));
            end
            subCarb = sum(subCarb);

            for j = 1:size(prdLM,1)
                prdCarb = carbonCount.carbonAtoms(prdLM,1)*...
                    abs(full(model.S(prdLM(j,1),lrgMlcRxns(i,1))));
            end
            prdCarb = sum(prdCarb);
        else
            if size(subLM,1) > 1
                for j = 1:size(subLM,1)
                    subCarb(j,1) = carbonCount.carbonAtoms(subLM(j,1),1)*...
                        abs(full(model.S(subLM(j,1),lrgMlcRxns(i,1))));
                end
                subCarb = sum(subCarb);
                prdCarb = 0;
            else
                subCarb = carbonCount.carbonAtoms(subLM(1,1),1)*...
                    abs(full(model.S(subLM(1,1),lrgMlcRxns(i,1))));
                prdCarb = 0;
            end 
        end
    elseif size(subLM,1) < size(prdLM,1)
        if size(subLM,1) > 1
            for j = 1:size(prdLM,1)
                prdCarb(j,1) = carbonCount.carbonAtoms(prdLM(j,1),1)*...
                    abs(full(model.S(prdLM(j,1),lrgMlcRxns(i,1))));
            end
            prdCarb = sum(prdCarb);

            for j = 1:size(subLM,1)
                subCarb = carbonCount.carbonAtoms(subLM(j,1),1)*...
                    abs(full(model.S(subLM(j,1),lrgMlcRxns(i,1))));
            end
            subCarb = sum(subCarb);
        else
            if size(prdLM,1) > 1
                for j = 1:size(prdLM,1)
                    prdCarb(j,1) = carbonCount.carbonAtoms(prdLM(j,1),1)*...
                        abs(full(model.S(prdLM(j,1),lrgMlcRxns(i,1))));
                end
                prdCarb = sum(subCarb);
                subCarb = 0;
            else
                prdCarb = carbonCount.carbonAtoms(prdLM(1,1),1)*...
                    abs(full(model.S(prdLM(1,1),lrgMlcRxns(i,1))));
                subCarb = 0;
            end 
        end
    elseif size(prdLM,1) == size(subLM,1)
        prdCarb = carbonCount.carbonAtoms(prdLM(1,1),1)*...
            abs(full(model.S(prdLM(1,1),lrgMlcRxns(i,1))));
        subCarb = carbonCount.carbonAtoms(subLM,1)*...
            abs(full(model.S(subLM(1,1),lrgMlcRxns(i,1))));
    end

    subCarbon(i,1) = subCarb;
    prdCarbon(i,1) = prdCarb;    
end

missMatches = find(subCarbon - prdCarbon);
lmRxnIDs = lrgMlcRxns;
lmRxnIDs(missMatches) = [];
subCarbon(missMatches) = [];
prdCarbon(missMatches) = [];

lrgMlcs.tiles = '*** Cyt/Qui/Heme Rxns ***';
lrgMlcs.time = date();
lrgMlcs.lmMetIDs = lrgIDs;
lrgMlcs.lmRxnIDs = lmRxnIDs;
lrgMlcs.subCarbon = subCarbon;
lrgMlcs.prdCarbon = prdCarbon;
if ~isempty(missMatches)
    lrgMlcs.missMatch = true;
    lrgMlcs.numberOfMissMatches = size(missMatches,1);
    lrgMlcs.mM_IDs = lrgMlcRxns(missMatches);
else
    lrgMlcs.missMatch = false;
end  
