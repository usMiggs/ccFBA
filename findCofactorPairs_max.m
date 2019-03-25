function cofactorPairs = findCofactorPairs_max(model,metPairs,E,e)
 
% pattern = {'nadh';'nadph';'atp';'gtp';'ctp';'utp';'ttp';'fadh2'};
pattern = {'nadh';'nadph';'atp';'gtp';'fadh2'};

if exist('model','var')
    if isempty(model)
        error('No model input')
    else
        fnames = fieldnames(model);
        if ~ismember('metFormulas',fnames)
            error('field metFormulas in cbmodel required')
        end
    end
else
    error('No model input')
end

 
if exist('metPairs','var')
    if isempty(metPairs)
        smallMets = findSmallMets(model);
        metPairs = countMetPairs_max(model,smallMets.smallMetIDs,true,'_');
    end
else
    smallMets = findSmallMets(model);
    metPairs = countMetPairs_max(model,smallMets.smallMetIDs,true,'_');
end

if ~ismember('E',fnames)
    [E,e] = constructElementalMatrix(model.metFormulas);
    E = full(E);
else
    if ~ismember('elements',fnames)
        [E,e] = constructElementalMatrix(model.metFormulas);
        E = full(E);
    else
        if issparse(model.E)
            E = full(model.E);
            e = model.elements;
        else
            E = model.E;
            e = model.elements;
        end
    end
end

%elements selected for analysis
h = find(ismember(e,'H'));
p = find(ismember(e,'P'));
s = find(ismember(e,'S'));
o = find(ismember(e,'O'));

eGroup = [h,p,o];

c = 0;
cP = [];
for i = 1:size(metPairs,1)
    metPair = metPairs{i,6};
    metPair = cell2mat(metPair);
    %if E(metPair(1,1),s) == 0
        if E(metPair(1,1),p) ~= 0
            vec = bsxfun(@minus,E(metPair(1,1),:),E(metPair(2,1),:));
            vec2 = find(vec);
            if size(vec2,2) < 3
                if ismember(vec2,eGroup)
                    c = c + 1;
                    cP(c,1) = i;
                end
            end
        end
    %end
end

if isempty(cP)
    cofactorPairs_tmp = [];
    warning('no cofactor Pairs found')
else
    cofactorPairs_tmp = metPairs(cP,:);
end

 
for i = 1:size(pattern,1)
    lM(:,i) = ~cellfun(@isempty,regexpi(cofactorPairs_tmp(:,1),pattern(i,1)));
end

lV = logical(sum(lM,2));

cofactorPairs = cofactorPairs_tmp(lV,:);
