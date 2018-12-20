function minmax = runMinMax_GF(model,start_rxn_index,end_rxn_index,verbose);
[num_mets num_rxns] = size(model.S);

if nargin < 3
    start_rxn_index=1;
    end_rxn_index=num_rxns;
    verbose=true;
end

progressbar
count=1;

for i=start_rxn_index:end_rxn_index %(num_rxns-3):num_rxns
   
   %if verbose && ~isfield(model,'CS_varNames')
        %fprintf('minmax for %s\t',model.rxns{i});
        %model = changeObjective(model,model.rxns{i});
   %elseif verbose && isfield(model,'CS_varNames')
        %fprintf('minmax for %s\t',model.CS_varNames{i});
        model.c=zeros(size(model.S,2),1);
        model.c(i)=1;
   %end 
     
   sol = optimizeCbModel(model,'max');
   
   max(count,1) = sol.f;
   
   sol = optimizeCbModel(model,'min');
   
   min(count,1) = sol.f;
   
   %if verbose
        %fprintf('min: %d\t max: %d\n',min(count,1),max(count,1));
   %end
   
   count=count+1;
   progressbar(i/size(model.rxns,1))
end

minmax(:,1)=min;
minmax(:,2)=max;