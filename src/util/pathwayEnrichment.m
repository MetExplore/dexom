function [sorted_SubSystEnrichment, SubSystEnrichment] = pathwayEnrichment(model,rxnsList)
% perform pathway enrichment (as implemented in MetExplore http://metexplore.fr/) from a list of reactions
% based on Fisher exact test, with Bonferroni correction (corrects for all performed tests, i.e., for all the susbsystems which include rxns from the input rxnsList)
% author: Nathalie POUPIN <nathalie.poupin at inrae.fr>

% find which model reactions are present in 'rxnsList'
InSampleRxns = zeros(length(model.rxns),1); % binary vector for model rxns, with 1 if reaction is in sample (in 'rxnsList') and 0 otherwise
for i = 1:length(rxnsList)
   InSampleRxns(strcmp(rxnsList{i},model.rxns),1)=1;
end

% find the number of reactions in each subsystem
[InSubsRxns,AllSubsystems] = countRxnsInSubsystems(model);
pvalue = zeros(size(InSubsRxns,2),1);

for i = 1:size(InSubsRxns,2)
    ContTable = crosstab(InSampleRxns,InSubsRxns(:,i));
    [~,p] = fishertest(ContTable,'Tail','right');
    pvalue(i) = p;
end

% Bonferroni correction
% -> corrects for all performed tests, i.e., for all the susbsystems which include rxns from the input rxnsList)
nbSubsWithRxns = sum(InSubsRxns'*InSampleRxns>0); 
    % NB: InSubsRxns'*InSampleRxns = nb. rxns from rxnsList are in each subSystem
pvalueCorr = pvalue*nbSubsWithRxns;
[~,ind] = sort(pvalueCorr);

nbRxnsPerSubsystem_inSample = (InSampleRxns'*InSubsRxns)';
nbRxnsPerSubsystem_inModel = sum(InSubsRxns,1)';
coverage = nbRxnsPerSubsystem_inSample./nbRxnsPerSubsystem_inModel;

FinalMatrix = [nbRxnsPerSubsystem_inModel nbRxnsPerSubsystem_inSample coverage pvalue pvalueCorr];
SubSystEnrichment = cat(2,AllSubsystems,num2cell(FinalMatrix));

sorted_FinalMatrix = FinalMatrix(ind,:);

sorted_subSystems = AllSubsystems(ind);
sorted_SubSystEnrichment = cat(2,sorted_subSystems,num2cell(sorted_FinalMatrix));
% sorted_SubSystEnrichment = cat(1,{'subsystems','total nb rxns','nb of mapped rxns','p-value','Bon. corrected p-value'},sorted_SubSystEnrichment);


function C=cleanCellArray(C)
i=1;
while i <= size(C,1)
    for j= 1:size(C,2)
        if iscell(C{i,j}) && size(C{i,j},1)<=1 && size(C{i,j},2)<=1 
            if ~isempty(C{i,j})
                C(i,j) = C{i,j};
            else
                 C{i,j}='';
            end
        end
        if isnumeric(C{i,j})
            c = C{i,j};
            c = (round(c*1000))/1000;
            C{i,j} = num2str(c);
        end
        if strcmp(C{i,j},'NaN')==1
            C{i,j}='';
        end
        if strcmp(C{i,j},'[]')==1 | strcmp(C{i,j},'''')==1
            C{i,j}='';
        end
        if strcmp(C{i,j},'unknown')==1
            C(i,:)='';
            i=i-1;
            break
        end
    end
    i=i+1;
end


function [RinS,uAllSubsystems] = countRxnsInSubsystems(model)

AllSubsystems = model.subSystems;
AllSubsystems = cleanCellArray(AllSubsystems);
uAllSubsystems = unique(AllSubsystems);

RinS = zeros(length(model.rxns),length(uAllSubsystems));
for j=1:length(uAllSubsystems)
    RinS(strcmp(AllSubsystems,uAllSubsystems{j}),j)=1;
end



