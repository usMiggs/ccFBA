function [model,isDrain] = resetModel(model,bnd,rev,zeroRxns)

%% Resets models bounds
%
%   *************************************************
%   *  CREATED by Maximilian Lularevic - 2017-05-24 *
%   *************************************************
%
% INPUTS
%   model   =   cobra model (needs model.rev,model.lb and model.ub fields)
%
%   bnd     =   any value to which you would like to set the
%               maximum/minimum bnds
%               default = 100
%
%   rev     =   a vector the size of model.rxns with 0 and 1 where 1 means
%               bidirectional and 0 unidirectional (such as model.rev).  if
%               bnds need to deviate from model.rev vector (tFBA)
%               default = model.rev
%
% OUTPUTS
%   model   =   model with reset lower and upper bounds
%
%   isDrain =   list of all drain rxns
%
%
%

if exist('bnd','var')
    if isempty(bnd)
        bnd = 100;
    end
else
    bnd = 100;
end

if exist('rev','var')
    if isempty(rev)
        rev = model.rev;
    end
    if size(rev) ~= size(model.rxns)
        error('rev vector and number of rxns do not match')
    end
else
    rev = model.rev;
end


% resetting model bnds
for i = 1:size(model.rxns)
    if rev(i,1) == 1
        model.lb(i,1) = -abs(bnd);
        model.ub(i,1) = abs(bnd);
    else
        model.lb(i,1) = 0;
        model.ub(i,1) = abs(bnd);
    end
end

% setting all uptake rxns to secretion only
isDrain = find(checkDrainRxns(model));
model.lb(isDrain) = 0;

if exist('zeroRxns','var')
    if ~isempty(zeroRxns)
        model.lb(zeroRxns) = 0;
        model.ub(zeroRxns) = 0;
    
        sol = optimizeCbModel(model);
        if ismember(sol.origStat,'INFEASIBLE')
            error('model does not solve with the zeroRxns')
        else
            strcat('RxnID:',num2str(zeroRxns),' ',' set to 0')
        end
    end
end


            