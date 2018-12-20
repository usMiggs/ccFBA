function model = setMinMax(model,minmax,tFBA)
%% This function is designed to one-line setting lb and ub according to the minmax
%
% INPUTS
%   model   =   cobra model (needs model.lb and model.ub fields)
%
%   minmax  =   minmax variable the size of model.rxns-by-2 where the first
%               column of the minmax matrix represents the model.lb and the
%               second column model.ub.
%
% OUTPUTS
%   model   =   cobra model with its bounds set to the minmax result input
%
%
%
%   2017-05-26 - Maximilian Lularevic
%

if nargin < 2
    error('not enough input arguements - cobra model and minmax matrix need to be provided')
end

if size(minmax,1) ~= size(model.rxns,1)
    error('dimension missmatch - make sure the minmax matrix has as many columns as the cobra model has rxns')
end

if ~isfield(model,'lb')
    error('cobra model structure for lb and ub not provided')
end

if exist('tFBA','var')
    if tFBA
        if isfield(model,'A')
            model = setTFBAminmax(model,minmax,'FBA')
            msgbox('netfluxes for var_lb and var_ub were set to input minmax in addition to regular lb and ub!')
        else
            model.lb = minmax(:,1);
            model.ub = minmax(:,2);
        end
    else
        model.lb = minmax(:,1);
        model.ub = minmax(:,2);
    end
else
    model.lb = minmax(:,1);
    model.ub = minmax(:,2);
end


    



