function isTrans = findTransports(model,rxns,delimiter)
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

% evaluate if specific reaction is looked for or all in model
if exist('rxns','var')
    if isempty(rxns)
        rxns = [1:size(model.rxns)]';
    else
        if iscellstr(rxns)
            rxns = find(ismember(model,rxns));
        elseif isstring(rxns)
            rxns = find(ismember(model,rxns));
        end
    end
else
    rxns = [1:size(model.rxns)]';
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
for i=1:length(rxns)
    reactantsID = find(model.S(:,rxns(i)) < 0);
    productsID = find(model.S(:,rxns(i)) > 0);
    IDs = [reactantsID;productsID];
    comps = {};
    
%     reactants = model.mets(find(model.S(:,rxns(i)) < 0));
%     products = model.mets(find(model.S(:,rxns(i)) > 0));

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
        if delimiter == '_'
            for j = 1:size(find(model.S(:,rxns(i))),1)
                met = model.mets{IDs(j,1)};
                comps{j,1} = met(end);
            end
            if size(unique(comps),1) > 1
                isTrans(i,1) = 1;
            else
                isTrans(i,1) = 0;
            end
        elseif delimiter =='['
            for j = 1:size(find(model.S(:,rxns(i))),1)
                met = model.mets{IDs(j,1)};
                comps{j,1} = substr(met,-2,1);
            end
            if size(unique(comps),1) > 1
                isTrans(i,1) = 1;
            else
                isTrans(i,1) = 0;
            end
        else
            error('delimiter does not match any predefined values. In order to prevent any further errors please add a metCompSymbol field to your model containing all metabolite compartment symbols.')
        end
    end
end

    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
