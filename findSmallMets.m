function [metInfo] = findSmallMets(model,NoFormulaStr,cut_off)
% Simple function to find inorganic metabolites from COBRA compatible
% models that have a metFormula entry (model.metFormulas).
%
% INPUTS

% evaluate whether NoFormulaStr was passed on (i.e. what is in 
% model.metFormula when no metFormula is present)
if exist('NoFormulaStr','var')
    if isempty(NoFormulaStr)
        NoFormulaStr = {'NA','None','NULL',''};
    end
else
    NoFormulaStr = {'NA','None','NULL',''};
end
    
% evaluate whether cut_off variable was passed on 
% cut_off value will ignore metabolites whose met string is longer 
% (i.e. has more characters)than the specified value
if exist('cut_off','var')
    if isempty(cut_off)
        cut_off = 7;
    end
else
    cut_off = 7;
end

% Make sure that the provided model contains a field with the chemical formulas
if ~isfield(model,'metFormulas')
    error('Model File Does not contain Metabolite Formulas or the variable is not named model.metFormulas')
end

% Find the indices of the metabolites whose chemical formulas does not
% contain 'C' followed by another capital letter or any number. For example
% indices of elements such as 'Ca', 'Cd', etc will be found:
InorganicMetID = find(cellfun(@isempty,regexp(model.metFormulas,'(?-i)C[A-Z_0-9]')));

% From these indices (InorganicMetID) exclude the indices that correspond
% to "NoFormulaStr".
NoFormulaStrID = find(ismember(model.metFormulas,NoFormulaStr));
InorganicMetID = setdiff(InorganicMetID,NoFormulaStrID);

% based on these indices (InorganicMetID) remove metIDs whose met string
% has more characters than the set cut_off value
for i = 1:size(InorganicMetID,1)
    logVec(i,1) = length(char(model.mets(InorganicMetID(i,1)))) < cut_off+1;
end
InorganicMetID = InorganicMetID(logVec);

% remove all lipids and glycos which are denoted as 'X' or at least contain
% a capital X in the metFormula (e.g. 'HX')
InorganicMetID=InorganicMetID(find(cellfun(@isempty,...
    regexp(model.metFormulas(InorganicMetID),'X'))),1);

% compile info in table
smallMetList = table(InorganicMetID,model.mets(InorganicMetID),...
    model.metFormulas(InorganicMetID),...
    'VariableNames',{'metID','metName','chemFormula'});

% Append all useful properties in structure
metInfo.title = '*** finding small mets ***';
metInfo.dateAndTime = datetime();
metInfo.cut_off = cut_off;
if isfield(model,'description')
    metInfo.model = model.description;
end
metInfo.mets = model.mets(InorganicMetID);
metInfo.metFormulas = model.metFormulas(InorganicMetID);
metInfo.smallMetIDs = InorganicMetID;
metInfo.table = smallMetList;



      