%% Oude Lohuis et al. 2023 Nat Neurosci
% Triple dissociation of auditory, visual, and motor processing in primary visual cortex
% MOL (C) 2023

% This script analyzes the cross prediction power of GLM on trial types of a different sort
% in the audiovisual change detection task

%%
startover

%% Load dataset:
% savedate                    = '03-01-23cross';
% folderpath                  = fullfile('E:','Data','Analysis','neuroGLM',savedate);
folderpath                  = fullfile('E:','Matlab','oudelohuis-et-al-2023-natneurosci','Fig3,4','GLMfits');
fileList                    = dir(fullfile(folderpath,'*.mat'));
fileList                    = {fileList(:).name};
fileList                    = fileList(~contains(fileList,'X'));
nFiles                      = length(fileList);

%Load the first session, then concatenate the rest
loadstruct          = load(fullfile(folderpath,fileList{1}));

output.x_sesid      = cell(nFiles,1);
output.y            = NaN(2000,size(loadstruct.output.y,2));
output.modelFits    = loadstruct.output.modelFits;
output.x_label      = loadstruct.output.x_label;
nTotalneurons       = 0;

sessionData         = struct();
trialData           = struct();
spikeData           = struct();

cvR2                = [];
cvR2_cat            = [];
cvR2_cross          = [];
cvR2_cross_cat      = [];

for iF = 1:nFiles
    fprintf('Loading and appending data from session #%d/%d\n',iF,nFiles)
    loadstruct          = load(fullfile(folderpath,fileList{iF}));

    sessionData         = AppendStruct(sessionData,loadstruct.sessionData);
    trialData           = AppendStruct(trialData,loadstruct.trialData);
    spikeData           = AppendStruct(spikeData,loadstruct.spikeData);
    
    nNeurons = length(loadstruct.spikeData.session_ID);
    idx = nTotalneurons+1:nTotalneurons+nNeurons;
    
    output.x_sesid(iF,1)= loadstruct.output.x_sesid;
    output.y(idx,:)     = loadstruct.output.y;
    if iF>1
        output.modelFits    = [output.modelFits; loadstruct.output.modelFits];
    end
    
    cvR2            = [cvR2; loadstruct.cvR2]; %#ok<*AGROW>
    cvR2_cat        = [cvR2_cat; loadstruct.cvR2_cat];
    cvR2_cross      = [cvR2_cross; loadstruct.cvR2_cross];
    cvR2_cross_cat  = [cvR2_cross_cat; loadstruct.cvR2_cross_cat];
    
    nTotalneurons       = nTotalneurons+nNeurons; 
    
end

sessionData.Experiment      = strrep(sessionData.Experiment,num2str(2),'');

%% Parameters:
params                      = loadstruct.params;
params                      = MOL_getColors_CHDET(params);
params.savedir              = 'E:\Documents\PhD\Figures\Project CHDET\Results - auV1\18GLM\';
params.Yvarlabels           = {'True' 'Full model' 'Trial' 'Visual' 'Audio' 'Motor'};

%% Remove animal 1012 with recordings in LM according to histology
spikeData.area(strcmp(spikeData.session_ID,'10122019041022') & strcmp(spikeData.area,'A1')) = deal({'LM'});
spikeData.area(strcmp(spikeData.session_ID,'10122019041126') & strcmp(spikeData.area,'A1')) = deal({'LM'});

%% Report dataset:
nSessions           = length(sessionData.session_ID);
nNeurons            = length(spikeData.session_ID);
nTrials             = length(trialData.session_ID);

output.y            = output.y(1:nNeurons,:);
fprintf('Dataset: %d sessions, %d trials, %d neurons, %d videos\n',nSessions,nTrials,nNeurons,length(sessionData.session_ID));

for iExp = 1:3
    fprintf('%s: %d sessions\n',params.ExperimentLabels{iExp},sum(strcmp(sessionData.Experiment,params.Experiments(iExp))));
    idxV1 = strcmp(spikeData.area,'V1') & ismember(spikeData.session_ID,sessionData.session_ID(strcmp(sessionData.Experiment,params.Experiments(iExp))));
    idxAC = strcmp(spikeData.area,'AC') & ismember(spikeData.session_ID,sessionData.session_ID(strcmp(sessionData.Experiment,params.Experiments(iExp))));
    fprintf('%d V1 and %d AC neurons \n',sum(idxV1),sum(idxAC));
    idx = ismember(trialData.session_ID,sessionData.session_ID(strcmp(sessionData.Experiment,params.Experiments(iExp))));
    fprintf('%d trials \n\n',sum(idx));
end

%% Multi-level statistics: %get var that indexes which mouse the session belongs to
G_mou           = cell(nNeurons,1);
uMice           = unique(sessionData.mousename);
for iMouse = 1:length(uMice)
    G_mou(ismember(spikeData.session_ID,sessionData.session_ID(strcmp(sessionData.mousename,uMice{iMouse})))) = uMice(iMouse);
end

%% Correct for negative R2: 
cvR2(cvR2<0)                       = NaN;
cvR2_cat(cvR2_cat<0)               = NaN;
cvR2_cross(cvR2_cross<0)           = NaN;
cvR2_cross_cat(cvR2_cross_cat<0)   = NaN;

%% Report explained variance in this dataset with conflict trials etc.
idx     = strcmp(spikeData.area,'V1');
fprintf('Variance explained across all V1 neurons (%2.3f, IQR %2.3f-%2.3f)\n',...
    nanmedian(cvR2(idx,1)),prctile(cvR2(idx,1),[25 75]));

%% FIGURES: 

%% Fig 6g: showing for each predictor the translative power:
idx_exp       = ismember(spikeData.session_ID,sessionData.session_ID(strcmp(sessionData.Experiment,'ChangeDetectionConflict')));
idx           = strcmp(spikeData.area,'V1') & idx_exp; %take only V1 and only MST animals

temp          = cvR2_cross_cat(idx,:,:,:);
temp          = squeeze(nanmean(temp,1));

params.trialTypes       = {'C' 'V' 'A' 'AV'};
params.labels_cats      = {'Raw' 'Full' 'Trial' 'Vis' 'Aud' 'Motor'};

cmap = parula(128);
cmap = [0.7 0.7 0.7; cmap]; %force 0 cvR2 to gray color

figure; set(gcf,'units','normalized','Position',[0.05 0.6 0.7 0.21],'color','w'); hold all;
for iC = 1:4
    subplot(1,4,iC); hold all;
    tmp     = temp(:,:,iC);
    imagesc(tmp)
    colormap(cmap)
    caxis([0 0.08])
    title(params.labels_cats{iC+2});
    set(gca,'XTick',1:4,'YTick',1:4,'XTickLabels',params.trialTypes,'YTickLabels',params.trialTypes)
    ylabel('Test')
    xlabel('Train')
    xlim([0.5 4.5])
    ylim([0.5 4.5])
    c = colorbar(); c.Ticks = [0 0.08];
end

export_fig(fullfile(params.savedir,sprintf('heatmap_cvR2_crosspred')),'-eps','-nocrop')

%% Statistics:
idx_exp         = ismember(spikeData.session_ID,sessionData.session_ID(strcmp(sessionData.Experiment,'ChangeDetectionConflict')));
idx_V1           = strcmp(spikeData.area,'V1') & idx_exp;

vis_within      = [cvR2_cross_cat(idx_V1,2,2,2); cvR2_cross_cat(idx_V1,4,4,2)];
vis_cross       = [cvR2_cross_cat(idx_V1,2,4,2); cvR2_cross_cat(idx_V1,4,2,2)];

X_var           = [ones(sum(idx_V1)*2,1); ones(sum(idx_V1)*2,1)*2];
G_mou2          = repmat(G_mou(idx_V1),4,1);

Y               = [vis_within; vis_cross];
tbl             = table(Y,X_var,G_mou2,'VariableNames',{'Y','Var','Mouse'}); %Create table for mixed model
lme             = fitlme(tbl,'Y~Var+(1|Mouse)'); %construct linear mixed effects model with fixed effect of temporal window (NO random intercept for mice as cohorts are tested)
stats           = dataset2table(anova(lme,'DFMethod','Satterthwaite')); %Perform ANOVA on model and output as matrix
fprintf('Visual: %1.3f vs %1.3f; F(%d,%2.0f)=%1.2f, p=%1.4f; \n',...
    nanmean(vis_within),nanmean(vis_cross),stats{2,3},stats{2,4},stats{2,2},stats{2,5})

writetable(tbl,'SourceData_Fig6f_GLM_CrossVis.xlsx')

aud_within        = [cvR2_cross_cat(idx_V1,3,3,3); cvR2_cross_cat(idx_V1,4,4,3)];
aud_cross         = [cvR2_cross_cat(idx_V1,3,4,3); cvR2_cross_cat(idx_V1,4,3,3)];

Y               = [aud_within; aud_cross];
tbl             = table(Y,X_var,G_mou2,'VariableNames',{'Y','Var','Mouse'}); %Create table for mixed model
lme             = fitlme(tbl,'Y~Var+(1|Mouse)'); %construct linear mixed effects model with fixed effect of temporal window (NO random intercept for mice as cohorts are tested)
stats           = dataset2table(anova(lme,'DFMethod','Satterthwaite')); %Perform ANOVA on model and output as matrix
fprintf('Auditory: %1.3f vs %1.3f; F(%d,%2.0f)=%1.2f, p=%1.2f; \n',...
    nanmean(aud_within),nanmean(aud_cross),stats{2,3},stats{2,4},stats{2,2},stats{2,5})

writetable(tbl,'SourceData_Fig6f_GLM_CrossAu.xlsx')

motor_within        = [cvR2_cross_cat(idx_V1,2,2,4); cvR2_cross_cat(idx_V1,3,3,4); cvR2_cross_cat(idx_V1,4,4,4)];
motor_cross         = [cvR2_cross_cat(idx_V1,2,4,4); cvR2_cross_cat(idx_V1,3,4,4); cvR2_cross_cat(idx_V1,4,3,4); cvR2_cross_cat(idx_V1,4,2,4)];

X_var           = [ones(sum(idx_V1)*3,1); ones(sum(idx_V1)*4,1)*2];
G_mou2          = repmat(G_mou(idx_V1),7,1);

Y               = [motor_within; motor_cross];
tbl             = table(Y,X_var,G_mou2,'VariableNames',{'Y','Var','Mouse'}); %Create table for mixed model
lme             = fitlme(tbl,'Y~Var+(1|Mouse)'); %construct linear mixed effects model with fixed effect of temporal window (NO random intercept for mice as cohorts are tested)
stats           = dataset2table(anova(lme,'DFMethod','Satterthwaite')); %Perform ANOVA on model and output as matrix
fprintf('Motor: %1.3f vs %1.3f; F(%d,%2.0f)=%1.2f, p=%1.2f; \n',...
    nanmean(motor_within),nanmean(motor_cross),stats{2,3},stats{2,4},stats{2,2},stats{2,5})

writetable(tbl,'SourceData_Fig6f_GLM_CrossMotor.xlsx')
