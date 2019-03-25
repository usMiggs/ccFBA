function cStruc = findMetCarbon(model)
%% This script is designed to determine the number of carbon atoms per metabolite
%   This function goes through all metabolites of a cobra model, and generates
%   a list with the number of carbon atoms of that metabolite
%
%               *****************************************
%               *   created by - Maximilian Lularevic
%               *              2017-09-30               *
%               *****************************************
%
% INPUTS
%   model   =   cobra nodel (needs to contain model.metFormulas field)
%
% OUTPUT
%   cStruc  =   structure containing the following fields
%               
%
%   Version History:
%       Date:           Comment:                        Initials:
%       2017-11-17      Updated function to be more     MLU  
%                       accurate in carbon count 
%                       (COFULLR2COFULLR1 etc.)
%
%       2019-03-22      Updated function to catch       MLU
%                       metFormulas such as
%                       'C2H4NO2RC2H2NOR' which 
%                       include two C2s and therfore
%                       have 4 carbons, code would
%                       fail as this was not 
%                       encountered before (happened
%                       with iMM904 model)
%
%

% chekc if metFormula is a field in model structure
if ~isfield(model,'metFormulas')
    error('model structure needs metFormula field')
end

% count the number of mets and stroe model.mets in mets variable
[met_num,~] = size(model.S);
mets = model.mets;

% remove the compartment information from mets (i.e. nadh_c -> nadh) and
% stroe new mets in metList
for i = 1:met_num
    met_tmp = regexp(mets(i),'_[a-z]{1,2}$','split');
    metList{i,1} = met_tmp{1,1}{1,1};
end

% find all fields with no metFormula in it
NoFormulaStr = {'NA','None','NULL',''};
NoFormulaStrID = find(ismember(model.metFormulas,NoFormulaStr));

% check for unique metabolites. before processing mets into metList
% metabolites such as nadh_c and nadh_m existed. since processing there are
% multiples in that list (i.e. nadh and nadh etc.)
% uList(ic) = metList; hence ic can be used to reconstruct the carbon
% vector
[uList,ia,ic] = unique(metList);

% extract metFormul;as into variable metForm
metForm = model.metFormulas(ia);

% loop through all unique metabolites (ia = number of unique mets)
for i = 1:size(ia,1)
    % checks for all capital C's in metForm folowed by a digit
    carbon_tmp = regexp(metForm(i,1),'(?<=C)\d+','match');
    % this might fail as there are metabolites such as CH3, where capital C
    % is followed by a letter
    
    % edited for the case of a met folula such as 'C2H4NO2RC2H2NOR' whihc
    % includes two C2 atoms as one belongs to the R group
    if size(carbon_tmp{1,1},2) > 1
        c_tmp = carbon_tmp;
        s_c_tmp = sum(str2double(carbon_tmp{1,1}));
        clear carbon_tmp
        carbon_tmp{1,1} = int2str(s_c_tmp);
    end
    if isempty(carbon_tmp{1,1})
        carbon_tmp{1,1} = 0;
    end
        % this regexp function checks if there is a capital C folowed by
        % any capital letter or number. This is to prevent fales positives
        % such as NaCl where capital C is followed by lowercase letter l
    isCarbon = regexpi(metForm(i,1),'(?-i)C[A-Z]');
    if isempty(isCarbon{1,1})
        carbon2 = 0;
    else
        carbon2 = size(isCarbon{1,1},2);
    end
    carbon_tmp = str2double(carbon_tmp{1,1});
    if isnan(carbon_tmp)
        carbon_tmp = 0;
    end
    carbon(i,1) = carbon_tmp + carbon2;
    %else
        %carbon(i,1) = str2double(carbon_tmp{1,1});
    %end
end

% metabolite name met formulas and number of carbon atoms in table
carbonList = table([1:met_num]',mets,metForm(ic),carbon(ic),...
    'VariableNames',{'metID','metName','chemFormula','carbons'});

% constructing output structure cStruc
cStruc.title = '*** Met Carbon Count ***';
cStruc.dateAndTime = datetime();
if isfield(model,'description')
    cStruc.model = model.description;
end
cStruc.metName = mets;
cStruc.metForm = metForm(ic);
cStruc.carbonAtoms = carbon(ic);
cStruc.table = carbonList;
















