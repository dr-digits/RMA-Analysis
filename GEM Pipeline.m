%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Installation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cd C:\Users\User\Documents\MATLAB\Human-GEM\code
HumanGEMInstaller.install
savepath

cd("C:\Users\User\Documents\MATLAB\RAVEN\installation")
checkInstallation

cd C:\Users\User\Documents\MATLAB\cobratoolbox
initCobraToolbox

changeCobraSolver("gurobi")
setRavenSolver("gurobi")

cd  "C:\Users\User\Documents\MATLAB\Human-GEM"
load("Human-GEM.mat");  % loads model as structure named "ihuman"
model = ravenCobraWrapper(ihuman);

ihuman = addBoundaryMets(ihuman); 
essentialTasks = parseTaskList("C:\Users\User\Documents\MATLAB\Human-GEM\data\metabolicTasks\metabolicTasks_Essential.txt");

checkTasks(ihuman, [], true, false, false, essentialTasks);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Generating context-specific model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cd C:\Users\User\Documents\MATLAB
TCGA_tpm  = readtable("Models_Universal.txt");
TCGA_tpm.Properties.VariableNames = ["ensemble_id","Sarcopenic_Cohort1","Healthy_Cohort1","Sarcopenic_Cohort2","Healthy_Cohort2","Sarcopenic_Universal","Healthy_Universal"];

[~, n] = size(TCGA_tpm);
numSamp = n-1;
Models_Universal = cell(numSamp, 1);

for i = 1:numSamp
data_struct.genes = cellstr(TCGA_tpm{:, 1});
data_struct.tissues =  TCGA_tpm.Properties.VariableNames(i+1);
data_struct.levels = TCGA_tpm{:, i+1};
data_struct.threshold = 1;
Models_Universal{i} = getINITModel2(ihuman, data_struct.tissues{1}, [],[], data_struct, [], true, [], true, true, essentialTasks, [],[]);
Models_Universal{i}.id = data_struct.tissues{1};
end


save("Models_Universal.mat", "Models_Universal")

checkTasks(Models_Universal{1, 1}  , [], true, false, false, essentialTasks);
checkTasks(Models_Universal{2, 1}  , [], true, false, false, essentialTasks);
checkTasks(Models_Universal{3, 1}  , [], true, false, false, essentialTasks);
checkTasks(Models_Universal{4, 1}  , [], true, false, false, essentialTasks);
checkTasks(Models_Universal{5, 1}  , [], true, false, false, essentialTasks);
checkTasks(Models_Universal{6, 1}  , [], true, false, false, essentialTasks);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% For structuraly model comparison
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
prepData = prepHumanModelForftINIT(ihuman, false, "C:\Users\User\Documents\MATLAB\Human-GEM\data\metabolicTasks\metabolicTasks_Essential.txt", "C:\Users\User\Documents\MATLAB\Human-GEM\model\reactions.tsv");
save("prepData.mat", "prepData")

baseModel = prepData.refModel;
compMat = false(length(baseModel.rxns), length(Models_Universal));

for i = 1:size(compMat,2)
    compMat(:,i) = ismember(baseModel.rxns,Models_Universal{i}.rxns);
end

rng(1);
proj_coords = tsne(double(compMat"), "Distance", "hamming", "NumDimensions", 2, "Exaggeration", 6, "Perplexity", 6);

model_ids = ["Sarcopenic_Cohort1","Healthy_Cohort1","Sarcopenic_Cohort2","Healthy_Cohort2","Sarcopenic_Universal","Healthy_Universal"];
models = Models_Universal;
model_ids = arrayfun(@(i) models{i}.id, (1:numel(models))", "UniformOutput", false);

res = compareMultipleModels(models);
clustergram(res.structComp, "Symmetric", false, "Colormap", "bone", "RowLabels", res.modelIDs, "ColumnLabels", res.modelIDs);

rxn2Dmap = tsne(res.reactions.matrix", "Distance", "hamming", "NumDimensions", 2, "Perplexity", 6);

scatter(rxn2Dmap(:,1), rxn2Dmap(:,2));
hold on
text(rxn2Dmap(:,1), rxn2Dmap(:,2), res.modelIDs);

useModels = res.modelIDs;
keep = ismember(res.modelIDs, useModels);
subMat = res.subsystems.matrix(:, keep);

subCoverage = (subMat - mean(subMat, 2)) .\ mean(subMat, 2) * 100;

inclSub = any(abs(subCoverage) > 25, 2);
subNames = res.subsystems.ID(inclSub);

cg = clustergram(subCoverage(inclSub,:), "Colormap", redbluecmap, "DisplayRange", 100, "rowLabels", subNames, "columnLabels", useModels, "ShowDendrogram", "OFF");

taskFileName = "C:\Users\User\Documents\MATLAB\Human-GEM\data\metabolicTasks\metabolicTasks_Full.txt";

res_func = compareMultipleModels(models(keep), false, false, [], true, taskFileName);

isDiff = ~all(res_func.funcComp.matrix == 0, 2) & ~all(res_func.funcComp.matrix == 1, 2); 
diffTasks = res_func.funcComp.tasks(isDiff);
spy(res_func.funcComp.matrix(isDiff,:), 30);

set(gca, "XTick", 1:numel(useModels), "XTickLabel", useModels, "XTickLabelRotation", 90, ...
    "YTick", 1:numel(diffTasks), "YTickLabel", diffTasks, "YAxisLocation", "right");
xlabel(gca, "");

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Run reporter metabolite analysis with DE genes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
modelS = simplifyModel(Models_Universal{3,1});

DEG_Res = readtable("DEG_Cohort2.csv");

genes1=DEG_Res.ensgene;
lgFC1=DEG_Res.log2FoldChange;
pvalues1=DEG_Res.pvalue;

outFileName1 = "reporterMetabolites_1.txt";
repMets1 = reporterMetabolites(modelS,genes1,pvalues1,true,outFileName1,lgFC1);

all_1=repMets1(1,:);
only_up_1=repMets1(2,:);
only_down_1=repMets1(3,:);
colnames= {"mets" "metNames" "metZScores" "metPValues" "metNGenes" "meanZ" "stdZ"};
all_1_T = table(all_1.mets,all_1.metNames,all_1.metZScores,all_1.metPValues,all_1.metNGenes,all_1.meanZ, all_1.stdZ, "VariableNames", colnames);
only_up_1_T = table(only_up_1.mets,only_up_1.metNames,only_up_1.metZScores,only_up_1.metPValues,only_up_1.metNGenes,only_up_1.meanZ, only_up_1.stdZ, "VariableNames", colnames);
only_down_1_T = table(only_down_1.mets,only_down_1.metNames,only_down_1.metZScores,only_down_1.metPValues,only_down_1.metNGenes,only_down_1.meanZ, only_down_1.stdZ, "VariableNames", colnames);

writetable(only_up_1_T,"only_up_BIRA.xlsx");
writetable(only_down_1_T,"only_down_BIRA.xlsx");

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Merging GEMs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function mergedModel = mergeConditionSpecificModels(model1,model2)
% mergedModel = mergeConditionSpecificModels(model1,model2)
% model1 will overlap model2 in all cases except for gene associations
% where the union is used
mergedModel.description = 'Merged model';
mergedModel.id = 'Universal Healthy Model';
mergedModel.annotation = model1.annotation;
mergedModel.rxns = union(model1.rxns,model2.rxns);
mergedModel.genes = union(model1.genes,model2.genes);
mergedModel.mets = union(model1.mets,model2.mets);
mergedModel.S = sparse(length(mergedModel.mets),length(mergedModel.rxns));
[A1,B1] = ismember(mergedModel.mets,model2.mets);
[A2,B2] = ismember(mergedModel.rxns,model2.rxns);
mergedModel.S(A1,A2) = model2.S(B1(A1),B2(A2));
[A1,B1] = ismember(mergedModel.mets,model1.mets);
[A2,B2] = ismember(mergedModel.rxns,model1.rxns);
mergedModel.S(A1,A2) = model1.S(B1(A1),B2(A2));
rxnGeneMat1 = sparse(length(mergedModel.rxns),length(mergedModel.genes));
rxnGeneMat2 = rxnGeneMat1;
[A1,B1] = ismember(mergedModel.rxns,model2.rxns);
[A2,B2] = ismember(mergedModel.genes,model2.genes);
rxnGeneMat1(A1,A2) = model2.rxnGeneMat(B1(A1),B2(A2));
[A1,B1] = ismember(mergedModel.rxns,model1.rxns);
[A2,B2] = ismember(mergedModel.genes,model1.genes);
rxnGeneMat2(A1,A2) = model1.rxnGeneMat(B1(A1),B2(A2));
mergedModel.rxnGeneMat = logical(rxnGeneMat1+rxnGeneMat2);
[A1,B1] = ismember(mergedModel.rxns,model1.rxns);
[A2,B2] = ismember(mergedModel.rxns,model2.rxns);
mergedModel.ub = zeros(length(mergedModel.rxns),1);
mergedModel.ub(A2) = model2.ub(B2(A2));
mergedModel.ub(A1) = model1.ub(B1(A1));
mergedModel.lb = zeros(length(mergedModel.rxns),1);
mergedModel.lb(A2) = model2.lb(B2(A2));
mergedModel.lb(A1) = model1.lb(B1(A1));

mergedModel.grRules = cell(length(mergedModel.rxns),1);
for i = 1:length(mergedModel.grRules)
    mergedModel.grRules{i} = strjoin(mergedModel.genes(mergedModel.rxnGeneMat(i,:)),' or ');
end

% mergedModel.rxnComps = zeros(length(mergedModel.rxns),1);
% mergedModel.rxnComps(A2) = model2.rxnComps(B2(A2));
% mergedModel.rxnComps(A1) = model1.rxnComps(B1(A1));
mergedModel.rxnNames = cell(length(mergedModel.rxns),1);
mergedModel.rxnNames(A2) = model2.rxnNames(B2(A2));
mergedModel.rxnNames(A1) = model1.rxnNames(B1(A1));
mergedModel.subSystems = cell(length(mergedModel.rxns),1);
mergedModel.subSystems(A2) = model2.subSystems(B2(A2));
mergedModel.subSystems(A1) = model1.subSystems(B1(A1));
mergedModel.eccodes = cell(length(mergedModel.rxns),1);
mergedModel.eccodes(A2) = model2.eccodes(B2(A2));
mergedModel.eccodes(A1) = model1.eccodes(B1(A1));
mergedModel.c = zeros(length(mergedModel.rxns),1);
mergedModel.c(A2) = model2.c(B2(A2));
mergedModel.c(A1) = model1.c(B1(A1));
mergedModel.rev = zeros(length(mergedModel.rxns),1);
mergedModel.rev(A2) = model2.rev(B2(A2));
mergedModel.rev(A1) = model1.rev(B1(A1));
[A1,B1] = ismember(mergedModel.mets,model1.mets);
[A2,B2] = ismember(mergedModel.mets,model2.mets);
if isfield(model1,'unconstrained')
    mergedModel.unconstrained = zeros(length(mergedModel.mets),1);
    mergedModel.unconstrained(A2) = model2.unconstrained(B2(A2));
    mergedModel.unconstrained(A1) = model1.unconstrained(B1(A1));
end
mergedModel.b = zeros(length(mergedModel.mets),1);
mergedModel.b(A2) = model2.b(B2(A2));
mergedModel.b(A1) = model1.b(B1(A1));
mergedModel.metNames = cell(length(mergedModel.mets),1);
mergedModel.metNames(A2) = model2.metNames(B2(A2));
mergedModel.metNames(A1) = model1.metNames(B1(A1));
mergedModel.metComps = zeros(length(mergedModel.mets),1);
mergedModel.metComps(A2) = model2.metComps(B2(A2));
mergedModel.metComps(A1) = model1.metComps(B1(A1));
mergedModel.metFormulas = cell(length(mergedModel.mets),1);
mergedModel.metFormulas(A2) = model2.metFormulas(B2(A2));
mergedModel.metFormulas(A1) = model1.metFormulas(B1(A1));
%mergedModel.metMiriams = cell(length(mergedModel.mets),1);
%mergedModel.metMiriams(A2) = model2.metMiriams(B2(A2));
%mergedModel.metMiriams(A1) = model1.metMiriams(B1(A1));
mergedModel.comps = model1.comps;
mergedModel.compNames = model1.compNames;
% mergedModel.compOutside = model1.compOutside;
% [A1,B1] = ismember(mergedModel.genes,model1.genes);
% [A2,B2] = ismember(mergedModel.genes,model2.genes);
% mergedModel.geneComps = zeros(length(mergedModel.genes),1);
% mergedModel.geneComps(A2) = model2.geneComps(B2(A2));
% mergedModel.geneComps(A1) = model1.geneComps(B1(A1));


end