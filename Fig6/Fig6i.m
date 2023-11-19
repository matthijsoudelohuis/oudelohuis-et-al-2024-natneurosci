%% Load exported orientation decoding scores from python: 
datadir         = 'E:\Matlab\oudelohuis-et-al-2024-natneurosci\Fig6\';
filename        = fullfile(datadir,'SourceData_Fig6i_SVM_coefficients.csv');
tbl             = readtable(filename);

%% Correlation audiovisual and visual trials orientation decoding weights:
nNeurons = size(tbl,1);
for iNeuron = 1:nNeurons
    g = sprintf('%d',tbl.unit_id(iNeuron,1));
    temp{iNeuron,1} = g(1:4);
end

tbl.animal_id = temp;
lme             = fitlme(tbl,'svm_coef_conf~svm_coef_vis+(1|animal_id)'); %construct linear mixed effects model with fixed effect of temporal window and random intercept for different mice
stats           = dataset2table(anova(lme,'DFMethod','Satterthwaite')); %Perform ANOVA on model and output as matrix
fprintf('Correlation SVM coefficients V and AV:\n')
fprintf('R=%1.2f,F(%d,%2.0f)=%1.2f, p=%1.3e; \n',sqrt(lme.Rsquared.Ordinary),stats{2,3},stats{2,4},stats{2,2},stats{2,5})
