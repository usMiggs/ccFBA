function all_met_pairs=findMPs_max(model,smallMetIDs,flag,type)
% This function goes through all reactions of a cobra model, and generates
% a list with all the possible metabolite pairs of reactants and products
% that participate in the same reaction.
%
% INPUT
%   model           =   cobra model
%   
%   smallMetIDs     =   this is a double containing a predefined list of rxnIDs
%                       of small metabolites. These will be excluded from the search of
%                       metabolite pairs.
%   
%   flag            =   metPair independent of compartment?(true or flase)
%                       default = false
%
% OUTPUT
% - all_met_pairs: This is a n-by-5 cell, where n is the number of met
%   pairs without repetition, found in the model. The 5 columns correspond
%   to:
%   (1) names of the metabolite pair separated by a colon punctuation
%   (2) number of occurence of this metabolite pair
%   (3) Formula of substrate
%   (4) Formula of product
%   (5) cell with rxnIDs corresponding to metPairs
%
%
%   ***** Modified by MaximilianLularevic - 2017-09-30 *****
%
%

if exist('flag','var')
    if isempty(flag)
        flag = false;
    end
else
    flag = false;
end

if exist('type','var')
    if isempty(type)
        type = '[';
        warning('type was set to default value = "["')
    end
else
    type = '[';
    warning('type was set to default value = "["')
end

% Get the number of reactions and metabolites in the model
[~,num_rxns]=size(model.S);

% Initialize variable
all_met_pairs={};

% Find all small metabolites according to smallMetIDs variable in the model
% and set their entire corresponding rows of the stoichiometric matrix to zero.
model.S(smallMetIDs,:)=0;

% Go through each reaction and find:
for i=1:num_rxns
    if isempty(find(model.S(:,i)))
        continue
    end 
    % Substrate/Product met names
    subIDs = find(model.S(:,i) < 0);
    substrates=model.mets(subIDs);
    prodIDs = find(model.S(:,i) > 0);
    products=model.mets(prodIDs);

    if flag
        if type == '['
            for s = 1:size(substrates,1)
                met = substrates{s,1};
                met_tmp = regexp(met,'\[\w','split');
                substrates{s,1} = met_tmp{1,1};
            end

            for p = 1:size(products,1)
                met = products{p,1};
                met_tmp = regexp(met,'\[\w','split');
                products{p,1} = met_tmp{1,1};
            end
        elseif type == '_'
            for s = 1:size(substrates,1)
                met = substrates{s,1};
                met_tmp = regexp(met,'\_\w','split');
                substrates{s,1} = met_tmp{1,1};
            end

            for p = 1:size(products,1)
                met = products{p,1};
                met_tmp = regexp(met,'\_\w','split');
                products{p,1} = met_tmp{1,1};
            end
        else
            warning('type varibale was not recognized and and therfore mets are not processed. Not processed means the compartment information is not cleaced off and mets are grouped by compartment.')
        end
    end
    
    % Substrate/Product Formulas
    substrates_Form=model.metFormulas(model.S(:,i) < 0);
    products_Form=model.metFormulas(model.S(:,i) > 0);

    % Check if this is a transport reaction of a metabolite between
    % compartments (if yes: isTrans=1)
    isTrans=findTransports(model,i,type);
    
    % If it is not a transport reaction, AND if there are more than 1
    % participating metabolites in this reaction (e.g.: this is to exclude
    % drain reactions)
    if ~isTrans && (length(find(model.S(:,i))) > 1)
        % Go through each substrate
        for j=1:length(substrates)
            % Go through each product
            for k=1:length(products)
                % Form two versions of every substrate and product pair of
                % names: subname:prodname, and prodname:subname
                current_met_pair1=strcat(substrates{j},':',products{k});
                current_met_pair2=strcat(products{k},':',substrates{j});

                % - metFormula of the substrate and product
                current_metForm_pair1a=substrates_Form{j};
                current_metForm_pair1b=products_Form{k};
                current_metID1a = subIDs(j);
                current_metID1b = prodIDs(k);
                
                % Now check if this metabolite pair exists already, in any
                % direction: subname:prodname, or prodname:subname
                if size(all_met_pairs,1) > 0
                    check1=find(ismember(all_met_pairs(:,1),current_met_pair1));
                    check2=find(ismember(all_met_pairs(:,1),current_met_pair2));
                % if not
                else
                    check1=[];
                    check2=[];
                end
                
                % Only if the pair has not been added in the list so far
                % with any of the two possible name configurations
                % (i.e.: check1 & check2 are empty), then:
                if isempty(check1) && isempty(check2)                   
                    % Add a new entry in the list:
                    % 'name pair' 'number of occurences' 'Formula sub' 'Formula prod'
                    all_met_pairs{size(all_met_pairs,1)+1,1} = current_met_pair1;
                    all_met_pairs{size(all_met_pairs,1),2}   = 1;
                    all_met_pairs{size(all_met_pairs,1),3}   = current_metForm_pair1a;
                    all_met_pairs{size(all_met_pairs,1),4}   = current_metForm_pair1b;
                    all_met_pairs{size(all_met_pairs,1),5}{1,1}   = i;
                    all_met_pairs{size(all_met_pairs,1),6}{1,1}   = current_metID1a;
                    all_met_pairs{size(all_met_pairs,1),6}{2,1}   = current_metID1b;
                    
                % If the name pair exists in the first configuration subname:prodname
                elseif ~isempty(check1)                    
                    % Increase the occurence of this met-pair configuration by 1
                    IDrow = size(all_met_pairs{check1,5},1);
                    all_met_pairs{check1,2}=all_met_pairs{check1,2}+1;
                    all_met_pairs{check1,5}{IDrow+1,1} = i;
                    if ~ismember(cell2mat(all_met_pairs{check1,6}),current_metID1a)
                        x = size(all_met_pairs{check1,6},1);
                        all_met_pairs{check1,6}{x+1,1} = current_metID1a;
                        all_met_pairs{check1,6}{x+2,1} = current_metID1b;
                    end
                    
                % If the name pair exists in the second configuration prodname:subname
                elseif ~isempty(check2)                    
                    % Increase the occurence of this met-pair configuration by 1
                    IDrow = size(all_met_pairs{check2,5},1);
                    all_met_pairs{check2,2}=all_met_pairs{check2,2}+1;
                    all_met_pairs{check2,5}{IDrow+1,1} = i;
                    if ~ismember(cell2mat(all_met_pairs{check2,6}),current_metID1a)
                        x = size(all_met_pairs{check2,6},1);
                        all_met_pairs{check2,6}{x+1,1} = current_metID1a;
                        all_met_pairs{check2,6}{x+2,1} = current_metID1b;
                    end                    
                end
            end
        end
    end
end

