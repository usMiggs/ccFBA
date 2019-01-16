function [carbonConst] = carbonConstraints(model,bnds,condition,carbonCount,cofPairs,minmax,runFinalMinMax,limEX,totalCarbon,intBnds,skipRxns)
%% This script is designed to constrain the solution space based on the 
%   number of carbon atoms made available to the model by the exchange
%   reactions. This is done in order to avoid physiologically meaningless
%   futile cycles.
%
%   CREATED 2017-09-30 by Maximilian Lularevic
%
% INPUT
%   model       =   cobra model (must contain model.metFormulas)
%
%   bnds        =   bounds or constraints as cell or double. for double 
%                   first row corresponds to rxnIDs, second row to lower 
%                   and the third row to upper bnds
%
% OPTIONAL
%
%   condition   =   string containing information as to which condition is 
%                   run (e.g. lactate producing conditon)
%
%   carbonCount =   carbonCount Structure as received form function 
%                   findMetCarbon.m (needs to contain
%                   carbonCount.carbonAtoms field). If this is not passed
%                   into the function it will create said structure by
%                   running the findMetCarbon function.
%
%   cofParis    =   cofParis structure as received from function
%                   findCofactorPairs.m. If not passed into the function it
%                   will create said structure by running cofPairs
%                   function.
%
%   minmax      =   FVA for model under experimental constraints. (needs 
%                   to be processed with fixMinMax.m). If not passed into
%                   the function it will create said structure by running
%                   cofPairs function. 
%
%   runFinalMM  =   logical true or false. if ture a final MinMax (FVA)
%                   will be run
%
%   limEX       =   value to be used to constrain drain rxns. This is
%                   supposed to be a single carbon constraint (e.g. left
%                   over carbon can be set to CO2 rate which is a single
%                   carbon. If flux for CO2 (closing the carbon by mass balance)
%                   is 0.5 mmol/gDCW*h this would be the rate to feed into
%                   the limEX variable).
%
%   intBnds     =   bnds used to constrain all unconstraint reactions.
%                   Default value is 100.
%
% OUTPUT
%   carbonConst =   structure containing the following fields
%                       - Title
%                       - Date and Time of job
%                       - constraintsName (Name of bnds variable: e.g. 'LacProd')
%                       - model.description
%                       - total carbon consumed
%               *** optional ***
%                       - actual carbon used for constraining (only when
%                         passed into function)
%               *** optional ***
%                       - report
%                           - cell with the lenght of model.rxns which 
%                             states whther reaction was balanced or unbalanced  
%                           - number of balanced rxns
%                           - number of unbalanced rxns
%                           - number of set rxns (aka constraints)
%                           - number of rxns with no flux
%                       - bndAdjusted
%                           - rawData (logical 2 colum matrix)
%                           - number of lower bounds adjusted by cc
%                             algorithm
%                           - number of upper bnds adjusted by cc algorithm
%                       - carbonCount
%                           - dateAndTime
%                           - modelDescription
%                           - metName
%                           - metFormula
%                           - carbon atoms for each metabolite
%                           - table with all the above information compiled
%                       - cofactorPairs (cell with mets and met IDs and
%                         occurance of these mets with corresponding
%                         rxnIDs)
%                       - original minmax (either the one passed into the
%                         function or herein calculated one)
%                       - newMinMax passed on carbon balance
%                       - solSpaceReduction size of new solution space
%                         compared to old one (calculatons based on minmax
%                         range comparison)
%                       - solutionNew (optimizeCbModel output with new bnds)
%                       - EXrxnsAdjusted -> title
%                       - limEX (single carbon limit passed into fucntion
%                         to constrain the Exchange reactions output)
%                       - number of EX bnds adjusted
%                       - solSpaceReductionEX in % compared to original
%                       - EXspaceReduction (how much has the space shrunk
%                         just comparing old EX ranges with new ones
%                       - post New MinMax -> title
%                       - postNewMinMax (new minmax run with runMinMax_GF
%                         and cc bnds)
%                       - solSpacereductionPostCC
%                       - solutionPostNew


%% check if minmax was passed into function and, run minmax if it hasn't

if exist('intBnds','var')
    if isempty(intBnds)
        intBnds = 100;
    end
else
    intBnds = 100;
end

[~,num_rxns] = size(model.S);
model = resetModel(model,intBnds);
model = loadConstraints(model,bnds,'FBA');
if exist('minmax','var')
    if isempty(minmax)
        minmax = runMinMax_GF(model);
        minmax = fixMinMax(minmax);
    end
else
    minmax = runMinMax_GF(model);
    minmax = fixMinMax(minmax);
end

%% check if condition was passed into function
if exist('condition','var')
    if isempty(condition)
        condition = 'n/a';
    end
else
    condition = 'n/a';
end

% check if final minmax should be run
if exist('runFinalMinMax','var')
    if isempty(runFinalMinMax)
        runFinalMinMax = 1;
    end 
else
    runFinalMinMax = 1;
end

% check if there are any rxns to be skipped
if exist('skipRxns','var')
    if isempty(skipRxns)
        skipRxns = [];
    else
        [s,z] = size(skipRxns);
        if z > s
            skipRxns = skipRxns';
        end
    end
else
    skipRxns = [];
end

%% extract lower and upper bnds from bnds variable
if iscell(bnds)
    k = 0;
    [num_vars,~] = size(bnds);
    for i=1:num_vars
        varName = bnds{i,1};
        var_index = find(ismember(model.rxns,varName));
        if isempty(var_index)
            fprintf('%s not found\n',bnds{i,1});
        else
            k = k+1;
            varID(k,1) = var_index;
            varMin(k,1) = bnds{i,2};
            varMax(k,1) = bnds{i,3};
        end
    end
elseif isnumeric(bnds)
    varID = bnds(:,1);
    varMin = bnds(:,2);
    varMax = bnds(:,3);
end
varID_tmp = varID;
%% find corresponding metIDs from varIDs
% first remove all bnd rxn that are not true drains (e.g. biomass, IgG formation)
isDrain = find(checkDrainRxns(model));
noExRxns = find(~ismember(varID,isDrain));
% adjust the previously retrieved variables
varID(noExRxns) = [];
varMin(noExRxns) = [];
varMax(noExRxns) = [];
% get all drains without the set bnds
isDrainNew = isDrain(~ismember(isDrain,varID));

%% check if number of carbon atoms per metabolite was passed into function
if exist('carbonCount','var')
    if isempty(carbonCount)
        carbonCount = findMetCarbon(model);
    end
else
    carbonCount = findMetCarbon(model);
end
% extract the carbon vector (length of model.mets)
metCarb = carbonCount.carbonAtoms;

%% multiply the lower and upper bnds witht their respective number of carbon
% atoms
for i = 1:size(varID,1)
    metID = find(model.S(:,varID(i,1)));
    atomCount(i,1) = sum(carbonCount.carbonAtoms(metID,1));
end
varMinCarb = varMin.*atomCount;

% varMaxCarb = varMax.*metCarb(metID);
if exist('totalCarbon','var')
    if isempty(totalCarbon)
        total_carbon_cons = abs(sum(varMinCarb(varMinCarb<0)));
    elseif ~isnumeric(totalCarbon)
        h = warndlg('Warning the totalCarbon input is not numeric. The totalCarbon will be calculated based on the lower bnds of the input variables. If you wish to continue press OK or else kill the code and input a numeric value for totalCarbon.','Problem with input');
        total_carbon_cons = abs(sum(varMinCarb(varMinCarb<0)));
    else
        total_carbon_cons = totalCarbon;
        actual_total_carbon_cons = abs(sum(varMinCarb(varMinCarb<0)));
    end
else
    total_carbon_cons = abs(sum(varMinCarb(varMinCarb<0)));
end

%% check if cofactorPairs have been passed on and set cofactors to 0 in model.S
% check for cofPair variable and/or calculate new
if exist('cofPairs','var')
    if isempty(cofPairs)
        smallMets = findSmallMets(model);
        metPairs = countMetPairs_max(model,smallMets.smallMetIDs,true,'_');
        cofPairs = findCofactorPairs_max(model,metPairs);
    end
else
    smallMets = findSmallMets(model);
    metPairs = countMetPairs_max(model,smallMets.smallMetIDs,true,'_');
    cofPairs = findCofactorPairs_max(model,metPairs);
end

% initialize variable
cofMetIDs = [];
% retrieve all the metabolite IDs for cofactors
for i = 1:size(cofPairs,1)
    cof_tmp = cell2mat(cofPairs{i,6});
    cofMetIDs = [cofMetIDs;cof_tmp];
end

% set cofMets to 0 in model.S
model_tmp = model;
model_tmp.S(cofMetIDs,:) = 0;

% find all rxns containing coa mets on both sides of the rxn
coa = findCoAs(model_tmp);

% find all rxns containing cytochromes, quinones, and hemes
lrgMlc = findLrgMlcs(model,carbonCount);

%% calculating new bnds algorithm
minmax_new = minmax;
mmAdj = zeros(num_rxns,2);
for i = 1:num_rxns
    % evaluate if rxn is part of the measured constraints
    if ismember(i,varID_tmp)
        bericht{i,1} = 'set_bnd';
        continue
    elseif ismember(i,skipRxns)
        bericht{i,1} = 'skipped rxn';
        continue
    else
        % check if reaction carries flux under normal conditions
        if sum(abs(minmax(i,:))) == 0
            bericht{i,1} = 'no_flux';
            continue
        elseif ismember(i,isDrain)
            bericht{i,1} = 'balanced - drain';
            subIDs = find(model_tmp.S(:,i) < 0);
            prodIDs = find(model_tmp.S(:,i) > 0);
            
            % account for stoichiometric coefficients and add up carbon
            subCarbon = sum(abs(full(model_tmp.S(subIDs,i))).*metCarb(subIDs,1));
            prodCarbon = sum(abs(full(model_tmp.S(prodIDs,i))).*metCarb(prodIDs,1));
            
            if subCarbon == 0 && prodCarbon == 0
                continue
            elseif subCarbon == 0
                constraint = abs(total_carbon_cons/prodCarbon);
            else
                constraint = abs(total_carbon_cons/subCarbon);
            end
            if model_tmp.rev(i,1) == 1
                if abs(minmax(i,1)) > constraint
                    minmax_new(i,1) = -constraint;
                    mmAdj(i,1) = 1;
                end
                if abs(minmax(i,2)) > constraint
                    minmax_new(i,2) = constraint;
                    mmAdj(i,2) = 1;
                end
            else
                if abs(minmax(i,2)) > constraint
                    if constraint > abs(minmax(i,1))
                        minmax_new(i,2) = constraint;
                        mmAdj(i,2) = 1;
                    end
                end
            end            
        else
            % find substrates and products 
            subIDs = find(model_tmp.S(:,i) < 0);
            prodIDs = find(model_tmp.S(:,i) > 0);
            
            % account for stoichiometric coefficients and add up carbon
            subCarbon = sum(abs(full(model_tmp.S(subIDs,i))).*metCarb(subIDs,1));
            prodCarbon = sum(abs(full(model_tmp.S(prodIDs,i))).*metCarb(prodIDs,1));

            % subtract 21 carbons form total sub and prd carbon
            % CoA = C21H32N7O16P3S
            if ismember(i,coa.coaRxnIDs)
                if coa.missMatch
                    warning('Stoichiometry of CoA is not matching. Smaller number is subtracted from carbon balance. Check reaction reaction manually.')
                    if coa.prdCarbon(ismember(coa.coaRxnIDs,i)) > coa.subCarbon(ismember(coa.coaRxnIDs,i))
                        subCarbon = subCarbon - coa.subCarbon(ismember(coa.coaRxnIDs,i));
                        prodCarbon = prodCarbon - coa.subCarbon(ismember(coa.coaRxnIDs,i));
                    else
                        subCarbon = subCarbon - coa.prdCarbon(ismember(coa.coaRxnIDs,i));
                        prodCarbon = prodCarbon - coa.prdCarbon(ismember(coa.coaRxnIDs,i));
                    end
                else
                    subCarbon = subCarbon - coa.subCarbon(ismember(coa.coaRxnIDs,i));
                    prodCarbon = prodCarbon - coa.subCarbon(ismember(coa.coaRxnIDs,i));
                end
            end            
            
            if ismember(i,lrgMlc.lmRxnIDs)
                subCarbon = subCarbon - (lrgMlc.subCarbon(ismember(lrgMlc.lmRxnIDs,i)));
                prodCarbon = prodCarbon - lrgMlc.prdCarbon(ismember(lrgMlc.lmRxnIDs,i));
            end
            
%             
%           % add up carbons
%           subCarbon = sum(metCarb(subIDs));
%           prodCarbon = sum(metCarb(prodIDs));
            
            if subCarbon == 0 && prodCarbon == 0
                bericht{i,1} = 'balanced';
                continue
%             elseif subCarbon == 0
%                 disp(strcat('carbon created in rxn ',num2str(i)));
%             elseif prodCarbon == 0
%                 disp(strcat('carbon drained in rxn ',num2str(i)));
            end
            
            % check carbon balancing
            if subCarbon == prodCarbon
                bericht{i,1} = 'balanced';
                constraint = abs(total_carbon_cons/subCarbon);
                if model_tmp.rev(i,1) == 1
                    if abs(minmax(i,1)) > constraint
                        minmax_new(i,1) = -constraint;
                        mmAdj(i,1) = 1;
                    end
                    if abs(minmax(i,2)) > constraint
                       minmax_new(i,2) = constraint;
                       mmAdj(i,2) = 1;
                    end
                else
                    if abs(minmax(i,2)) > constraint
                        if constraint > abs(minmax(i,1))
                            minmax_new(i,2) = constraint;
                            mmAdj(i,2) = 1;
                        end
                    end
                end     
            else % if carbon is not balanced
                % X_tmp = regexp(model.metFormulas(497),'X')
                X_sub = [];
                X_prod = [];
                for f = 1:(size(subIDs,1))
                    X_sub(f,1) = count(model.metFormulas(subIDs(f,1)),'X');
%                     X_tmp = regexpi(model.metFormulas(subIDs(f,1)),'X');         
%                     if isnan(X_sub(f,1))
%                         X_sub(f,1) = 0;
%                     end
                end
                for f = 1:(size(prodIDs,1))
                    X_prod(f,1) = count(model.metFormulas(prodIDs(f,1)),'X');
                end
                X_log = find([X_sub;X_prod]);
                
                if find(ismember(cofMetIDs,find(model.S(:,i))))
                    if mod(sum([subCarbon;prodCarbon;metCarb(cofMetIDs...
                        (ismember(cofMetIDs,find(model.S(:,i)))))]),2) == 0 % if the sum of all carbon is even then it is balanced (not bulletproof but if off by only one carbon it works)
                        bericht{i,1} = 'balanced';
                        if subCarbon > prodCarbon % if balanced then its go with the smaller one (constraining less). This is because there seems to be only 1 cofactor in this equation (e.g. atp -> amp, this is not picked up as a cofactor pair but atp is a cofactor) so the carbon of atp/amp should not be considered when constraining
                            constraint = abs(total_carbon_cons/prodCarbon);
                        else
                            constraint = abs(total_carbon_cons/subCarbon);
                        end
                        if model_tmp.rev(i,1) == 1
                            if abs(minmax(i,1)) > constraint
                                minmax_new(i,1) = -constraint;
                                mmAdj(i,1) = 1;
                            end
                            if abs(minmax(i,2)) > constraint
                               minmax_new(i,2) = constraint;
                               mmAdj(i,2) = 1;
                            end
                        else
                            if abs(minmax(i,2)) > constraint
                                if constraint > abs(minmax(i,1))
                                    minmax_new(i,2) = constraint;
                                    mmAdj(i,2) = 1;
                                end
                            end
                        end
                    else
                        bericht{i,1} = 'unbalanced - CofPair';
                        if subCarbon < prodCarbon % if unbalanced then go with the smaller one (constraining less)
                            constraint = abs(total_carbon_cons/subCarbon);
                        else
                            constraint = abs(total_carbon_cons/prodCarbon);
                        end
                        if model_tmp.rev(i,1) == 1
                            if abs(minmax(i,1)) > constraint
                                minmax_new(i,1) = -constraint;
                                mmAdj(i,1) = 1;
                            end
                            if abs(minmax(i,2)) > constraint
                               minmax_new(i,2) = constraint;
                               mmAdj(i,2) = 1;
                            end
                        else
                            if abs(minmax(i,2)) > constraint
                                if constraint > abs(minmax(i,1))
                                    minmax_new(i,2) = constraint;
                                    mmAdj(i,2) = 1;
                                end
                            end
                        end
                    end
                else
                    bericht{i,1} = 'unbalanced';
                    if subCarbon < prodCarbon % if unbalanced then go with the smaller one (constraining less)
                        constraint = abs(total_carbon_cons/subCarbon);
                    else
                        constraint = abs(total_carbon_cons/prodCarbon);
                    end
                    if model_tmp.rev(i,1) == 1
                        if abs(minmax(i,1)) > constraint
                            minmax_new(i,1) = -constraint;
                            mmAdj(i,1) = 1;
                        end
                        if abs(minmax(i,2)) > constraint
                           minmax_new(i,2) = constraint;
                           mmAdj(i,2) = 1;
                        end
                    else
                        if abs(minmax(i,2)) > constraint
                            if constraint > abs(minmax(i,1))
                                minmax_new(i,2) = constraint;
                                mmAdj(i,2) = 1;
                            end
                        end
                    end
                end
            end
        end
    end
end

%% qunatify how much solution space was constraint compared to original minmax
mm_diff = minmax(:,2)-minmax(:,1);
mm_new_diff = minmax_new(:,2)-minmax_new(:,1);
diff = 100*(1-(sum(mm_new_diff)/sum(mm_diff)));
        
%% variable complier
unbRxns = sum(ismember(bericht,'unbalanced'))+sum(ismember(bericht,'unbalanced - CofPair'));
balRxns = sum(ismember(bericht,'balanced'))+sum(ismember(bericht,'balanced - X'))...
    +sum(ismember(bericht,'balanced - drain'))+sum(ismember(bericht,'unconstr - oxPhos'));
setRxns = sum(ismember(bericht,'set_bnd'));
nfRxns = sum(ismember(bericht,'no_flux'));

carbonConst.title = '*** carbon closure constraining ***';
carbonConst.dateAndTime = datetime();
carbonConst.constraintsName = condition;
if isfield(model,'description')
    carbonConst.model = model.description;
end
carbonConst.totalCarbonConsumed = total_carbon_cons;
if exist('actual_total_carbon_cons','var')
    if ~isempty(actual_total_carbon_cons)
        carbonConst.actualTotalCarbonConsumed = actual_total_carbon_cons;
    end
end
carbonConst.report.summary = bericht;
carbonConst.report.unbalanced = unbRxns;
carbonConst.report.balanced = balRxns;
carbonConst.report.setRxns = setRxns;
carbonConst.report.noFlux = nfRxns;
carbonConst.bndAjusted.rawData = mmAdj;
carbonConst.bndAjusted.numberOf_LB_adj = sum(mmAdj(:,1));
carbonConst.bndAjusted.numberOf_UB_adj = sum(mmAdj(:,2));
carbonConst.CoAevaluation = coa;
carbonConst.electrTransferMlcs = lrgMlc;
carbonConst.carbonCount = carbonCount;
carbonConst.cofactorPairs = cofPairs;
carbonConst.originalMinMax = minmax;
carbonConst.newMinMax = minmax_new;
carbonConst.solSpaceReduction = strcat(num2str(diff),'%');

%% test if model solves under new conditions
model.lb = minmax_new(:,1);
model.ub = minmax_new(:,2);
sol = optimizeCbModel(model,'max');

if ismember(cellstr(sol.origStat),{'INFEASIBLE'})
    error('Model does not solve. Think about relaxing the bounds')
else
    disp(strcat('model solves. sol.x =',' ',num2str(sol.f)))
end

carbonConst.solutionNew = sol;

%% constrain uptakes based on passed in limit (e.g. lost/closed carbon)

if exist('limEX','var')
    if ~isempty(limEX)
        carbonConst.EXrxnAdjust = '*** EX rxns adjustment ***';
        carbonConst.limEX = limEX;
        
        % get all the max values for the EX rxns not set
        varMaxN = minmax_new(isDrainNew,2);
        % get met indices
        metID_ind = find(model.S(:,isDrainNew));
        % get metIDs from Met indices
        [metID_EX,~] = ind2sub(size(model.S),metID_ind);
        % get rxnIDs where rxn actually has a flux and ignore any fluxes
        % that are set to 0
        tmpRxnID = isDrainNew(find(varMaxN));
        % how many rxns have been set to 0?
        numEXzero = size(isDrainNew,1)-size(tmpRxnID,1);
        % get metIDs from EX rxn IDs so carbon per met can be determined
        tmpMetID = metID_EX(find(varMaxN));
        % calculate new bnds
        limEXnew = limEX./metCarb(tmpMetID);
        % get rid of Inf
        limEXnew(~isfinite(limEXnew)) = 0;
        % find which reactions needf adjustment
        adjEX = tmpRxnID(limEXnew < varMaxN(ismember(isDrainNew,tmpRxnID)));
        %  check if there is any bnds to be adjusted
        if ~isempty(adjEX)
            % make calculations as to how much solution space shrinks
            % compared to original
            diffLim = 100*(1-(sum(limEXnew(ismember(tmpRxnID,adjEX)))/sum(minmax_new(adjEX,2))));
            limEXadj = limEXnew(ismember(tmpRxnID,adjEX));
            minmax_new(adjEX,2) = limEXnew(ismember(tmpRxnID,adjEX));
            
            mm_lim_diff = minmax_new(:,2)-minmax_new(:,1);
            diff_n = 100*(1-(sum(mm_lim_diff)/sum(mm_diff)));
            
            model.lb = minmax_new(:,1);
            model.ub = minmax_new(:,2);
            sol = optimizeCbModel(model,'max');

            if isempty(sol.f)
                warning('model does not solve. think about relaxing the bounds')
            else
                disp(strcat('model solves. sol.x =',' ',num2str(sol.f)))
            end
            
            carbonConst.numAdjusted = size(adjEX,1);
            carbonConst.solSpaceReductionEX  = strcat(num2str(diff_n),'%');
            carbonConst.EXspaceReduction = strcat(num2str(diffLim),'%');
        else
            warning('no adjustment of EXrxns was done as already constrained more than limEX')
            carbonConst.warning = 'no adjustment was done';
        end
    end
end

%% run minmax on new minmax bnds and see if solution space is reduced further
if runFinalMinMax
    model.lb = minmax_new(:,1);
    model.ub = minmax_new(:,2);
    minmax_post_new = runMinMax_GF(model);
    carbonConst.postNmm = '*** post new minmax ***';
    carbonConst.postNewMinMax = minmax_post_new;
    % qunatify how much solution space was constraint compared to original minmax
    mm_post_diff = minmax_post_new(:,2)-minmax_post_new(:,1);
    diff = 100*(1-(sum(mm_post_diff)/sum(mm_diff)));
    carbonConst.solSpaceReductionPostCC = strcat(num2str(diff),'%');

    % test if new minmax solves
    model.lb = minmax_post_new(:,1);
    model.ub = minmax_post_new(:,2);
    sol = optimizeCbModel(model,'max');

    if isempty(sol.f)
        warning('model does not solve. think about relaxing the bounds')
    else
        disp(strcat('model solves. sol.x =',' ',num2str(sol.f)))
    end

    carbonConst.solutionPostNew = sol;
end

