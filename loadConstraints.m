function model = loadConstraints(model,constraints,type,print)
%% Function will constrain model based on input matrix
%
%   **************************************************
%   *  MODIFIED by Maximilian Lularevic - 2017-06-21 *
%   **************************************************
%
% INPUT:
%   model       =   COBRA model
%   constraints =   bnds matrix -> either *cell array* or *double*
%   type        =   'FBA' or 'TFBA' or 'NF-TFBA'(type of constarint) -> default 'FBA'
%   print       =   true or false -> if true will print out each constraint
%                   that was set
%
% OUTPUT:
%   model       =   COBRA model with set constraints. only rxns specified
%                   in the constraints matrix are altered
%
%  +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%  + WARNING: Double only works for FBA not TFBA constraints +
%  +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%
% Additional info for constraints:
% IF CELL:
%   following fields MUST be present:
%       -> first column specifies the variable names as per the varNames (or model.rxns for 'FBA')
%       -> second column specifies the lower bound
%       -> third colum specifies the upper bound
%
% IF DOUBLE:
%   following fields MUST be present:
%       -> first column specifies rxn idices according to model.rxns 
%       -> second column specifies the lower bound
%       -> third colum specifies the upper bound

%% function code
if iscell(constraints)
    if nargin < 3
        type = 'FBA';
    end

    if nargin < 4
        print = false;
    end

    if strcmp(type,'TFBA')
        [num_vars,~] = size(constraints);

        for i=1:num_vars
            varName = constraints{i,1};
            var_index = find(ismember(model.varNames,varName));
            if isempty(var_index)
                fprintf('%s not found\n',constraints{i,1});
            else
                if print
                    fprintf('%s lb: %d => %d\tub: %d => %d\n',constraints{i,1},model.var_lb(var_index),constraints{i,2},model.var_ub(var_index),constraints{i,3});
                end

                model.var_lb(var_index) = constraints{i,2};
                model.var_ub(var_index) = constraints{i,3};
            end
        end
    elseif strcmp(type,'NF-TFBA')
        [num_vars,~] = size(constraints);
        
        for i=1:num_vars
            varName = strcat('NF_',constraints{i,1});
            var_index = find(ismember(model.varNames,varName));
            if isempty(var_index)
                fprintf('%s not found\n',constraints{i,1});
            else
                if print
                    fprintf('%s lb: %d => %d\tub: %d => %d\n',constraints{i,1},model.var_lb(var_index),constraints{i,2},model.var_ub(var_index),constraints{i,3});
                end

                model.var_lb(var_index) = constraints{i,2};
                model.var_ub(var_index) = constraints{i,3};
            end
        end
            
    elseif strcmp(type,'FBA')
        [num_vars,~] = size(constraints);

        for i=1:num_vars
            varName = constraints{i,1};
            var_index = find(ismember(model.rxns,varName));

            if isempty(var_index)
                fprintf('%s not found\n',constraints{i,1});
            else
                if print
                    fprintf('%s lb: %d => %d\tub: %d => %d\n',constraints{i,1},model.lb(var_index),constraints{i,2},model.ub(var_index),constraints{i,3}); 
                end

                model.lb(var_index) = constraints{i,2};
                model.ub(var_index) = constraints{i,3};
            end
        end
    end
elseif isnumeric(constraints)
    for i = 1:size(constraints,1)
        model = changeRxnBounds(model,model.rxns(constraints(i,1)),constraints(i,2),'l');
        model = changeRxnBounds(model,model.rxns(constraints(i,1)),constraints(i,3),'u');
    end
else
    error('Constraints matrix does not match any know type')
end

        