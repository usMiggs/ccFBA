function isTrans = findTransports(model,delimiter)
% checks if the reaction is a transport reaction
% a transport reaction is defined as having the same metabolite(s) in different
% compartments

%% check input variable is correct
% check if cobra model was passed into function
if exist('model','var')
    if isempty(model)
        error('please supply cobra model')
    end
else
    error('please supply cobra model')
end

% evaluate if delimeter has been passed into function
if exist('delimiter','var')
    if isempty(delimiter)
        delimiter = '_';
        warning('default delimiter "_" was chosen!')
    end
else
    delimiter = '_';
    warning('default delimiter "_" was chosen!')
end

%% find trasport reactions code
for i=1:length(model.rxns)
    reactantsID = find(model.S(:,i) < 0);
    productsID = find(model.S(:,i) > 0);
    
    reactants = model.mets(find(model.S(:,i) < 0));
    products = model.mets(find(model.S(:,i) > 0));

    if isfield(model,'metCompSymbol')
        reac_comps = model.metCompSymbol(reactantsID);
        prod_comps = model.metCompSymbol(productsID);
        comps = [reac_comps;prod_comps];
        if size(unique(comps),1) > 1
            isTrans(i,1) = 1;
        else
            isTrans(i,1) = 0;
        end
    else
        for j = 1:size(reactantsID,1)+size(productsID,1)
            temp = regexp(fliplr(model.mets(reactantsID(j,1))),delimiter,'split');
            comps{j,1} = temp(1);
        end
        if size(unique(comps),1) > 1
            isTrans(i,1) = 1;
        else
            isTrans(i,1) = 0;
        end
    end
end

    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
