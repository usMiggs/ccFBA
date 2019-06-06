function isDrain=checkDrainRxns(model)
%% Modified checkDrainRxns function
%
% modified by Maximilian Lularevic 2018-08-13
%
% The function was modified in order to find drain reactions that contain
% more than one metabolite to drain, i.e proton coupled drains which is the
% case for some models such as the mitoCore model
%

if isfield(model,'isDrain')
    model = rmfield(model,'isDrain');
end

isDrain=zeros(length(model.rxns),1);

for i = 1:size(model.rxns,1)
    substrates = find(model.S(:,i) < 0);
    products = find(model.S(:,i) > 0);
    
    if isempty(substrates) || isempty(products)
        isDrain(i,1) = 1;
    else
        isDrain(i,1) = 0;
    end
end
