
%%
startover

%% Load dataset:

savedate                    = '22-12-22confl';
folderpath                  = fullfile('E:','Data','Analysis','neuroGLM',savedate);
fileList                    = dir(fullfile(folderpath,'*.mat'));
fileList                    = {fileList(:).name};
fileList                    = fileList(~contains(fileList,'X'));
nFiles                      = length(fileList);

%Load the first session, then concatenate the rest
loadstruct          = load(fullfile(folderpath,fileList{1}));

% output.x            = NaN(nFiles,size(loadstruct.output.x,2),size(loadstruct.output.x,3));
output.x_sesid      = cell(nFiles,1);
output.y            = NaN(2000,size(loadstruct.output.y,2));
output.modelFits    = loadstruct.output.modelFits;
% output.shuffleFits  = loadstruct.output.shuffleFits;
output.x_label      = loadstruct.output.x_label;
nTotalneurons       = 0;

sessionData         = struct();
trialData           = struct();
spikeData           = struct();
% videoData           = struct();

cvR2        = NaN(2000,4);
cvR2_cat    = NaN(2000,4,4);
cvR2_full   = NaN(2000,1);


for iF = 1:nFiles
    fprintf('Loading and appending data from session #%d/%d\n',iF,nFiles)
    loadstruct          = load(fullfile(folderpath,fileList{iF}));

    sessionData         = AppendStruct(sessionData,loadstruct.sessionData);
    trialData           = AppendStruct(trialData,loadstruct.trialData);
    spikeData           = AppendStruct(spikeData,loadstruct.spikeData);
%     videoData           = AppendStruct(videoData,loadstruct.videoData);
    
    nNeurons = length(loadstruct.spikeData.session_ID);
    idx = nTotalneurons+1:nTotalneurons+nNeurons;
    
%     output.x(iF,:,:)    = loadstruct.output.x;
    output.x_sesid(iF,1)= loadstruct.output.x_sesid;
    output.y(idx,:)     = loadstruct.output.y;
    if iF>1
        output.modelFits    = [output.modelFits; loadstruct.output.modelFits];
%         output.shuffleFits  = [output.shuffleFits; loadstruct.output.shuffleFits];
    end
    
    cvR2(idx,:) = loadstruct.cvR2;
    cvR2_cat(idx,:,:) = loadstruct.cvR2_cat;
    cvR2_full(idx,:) = loadstruct.cvR2_full;
    
    nTotalneurons       = nTotalneurons+nNeurons;    
end

sessionData.Experiment      = strrep(sessionData.Experiment,num2str(2),'');

%% Parameters:
params                      = loadstruct.params;
params                      = MOL_getColors_CHDET(params);
params.savedir              = 'E:\Documents\PhD\Figures\Project CHDET\Results - auV1\18GLM\';
params.Yvarlabels           = {'True' 'Full model' 'Trial' 'Visual' 'Audio' 'Motor'};

%% Report dataset:
nSessions           = length(sessionData.session_ID);
nNeurons            = length(spikeData.session_ID);
nTrials             = length(trialData.session_ID);

output.y            = output.y(1:nNeurons,:);
fprintf('Dataset: %d sessions, %d trials, %d neurons, %d videos\n',nSessions,nTrials,nNeurons,length(sessionData.session_ID));

for iExp = 1:3
    fprintf('%s: %d\n',params.ExperimentLabels{iExp},sum(strcmp(sessionData.Experiment,params.Experiments(iExp))));
end
fprintf('%d V1 and %d AC neurons \n',sum(strcmp(spikeData.area,'V1')),sum(strcmp(spikeData.area,'A1')));

%%
params.lambdastring         = 'lambda_min';
% params.lambdastring         = 'lambda_1se';

%% Multi-level statistics: 
X_coh           = cell(nNeurons,1);
X_coh(ismember(spikeData.session_ID,sessionData.session_ID(ismember(sessionData.Experiment,params.Experiments{1}))))         = params.ExperimentLabels(1);
X_coh(ismember(spikeData.session_ID,sessionData.session_ID(ismember(sessionData.Experiment,params.Experiments{2}))))         = params.ExperimentLabels(2);
X_coh(ismember(spikeData.session_ID,sessionData.session_ID(ismember(sessionData.Experiment,params.Experiments{3}))))         = params.ExperimentLabels(3);

G_mou           = cell(nNeurons,1);
uMice           = unique(sessionData.mousename);
for iMouse = 1:length(uMice)
    G_mou(ismember(spikeData.session_ID,sessionData.session_ID(strcmp(sessionData.mousename,uMice{iMouse})))) = uMice(iMouse);
end

%%

% figure; 
% 
% nanmean(cvR2_full)
% nanmean(cvR2,1)
% nanmean(cvR2_cat,1)
% 
% (idx,:) = loadstruct.cvR2;
%     cvR2_cat(idx,:,:) = loadstruct.cvR2_cat;
%     cvR2_full(idx,:) = loadstruct.cvR2_full;



%% FIGURES: 

%% Show model prediction:
close all;
cell_IDs            = {};

%visually driven example cell:
% cell_IDs{end+1}     = '20442021042411333'; %nice visual example

%Auditory + motor driven example cell:
cell_IDs{end+1}     = '20122018081431146';
%Auditory no motor  example cell:
cell_IDs{end+1}     = '10122019041031329'; 
%Motor no aud example cell:
cell_IDs{end+1}     = '20122018081311186';

% cell_IDs{end+1}     = '20122018081311186'; %perhaps
% cell_IDs{end+1}     = '20122018081431167'; %au only
% cell_IDs{end+1}     = '20122018081431176'; %perpaps but low firing
% cell_IDs{end+1}     = '20302020011621092';
% cell_IDs{end+1}     = '20112018081011125'; %all motor
% cell_IDs{end+1}     = '20442021042811300';
% cell_IDs{end+1}     = '20442021042811376';
% cell_IDs{end+1}     = '20302020011621372';
% cell_IDs{end+1}     = '10092019030831134'; %Very visually driven, plus one condition au, no movement

params.exportfig    = 0;

lastsesid = [];

% for iNeuron = 1:5%params.nNeurons %loop over neurons
for iNeuron = find(ismember(spikeData.cell_ID,cell_IDs))'
    nTotalSpikeBins     = sum(~isnan(output.y(iNeuron,:)));
    nTrials             = nTotalSpikeBins / params.nTimebins;
%     idx_ses             = strcmp(sessionData.session_ID,spikeData.session_ID(iNeuron));
%     temptrialData       = MOL_getTempPerSes(sessionData.session_ID(idx_ses),trialData);
    temptrialData       = MOL_getTempPerSes(spikeData.session_ID(iNeuron),trialData);
    
    if ~(strcmp(spikeData.session_ID(iNeuron),lastsesid))
        load(fullfile(folderpath,sprintf('X_Full_GLM_%s.mat',spikeData.session_ID{iNeuron})))
%         X_full              = squeeze(x(1,1:nTotalSpikeBins,1:params.nPredictors));
        lastsesid           = spikeData.session_ID(iNeuron);
    end
    
%     X_full              = squeeze(output.x(idx_ses,1:nTotalSpikeBins,1:params.nPredictors));
    Y_full              = squeeze(output.y(iNeuron,1:nTotalSpikeBins))';
    
    Yh_full             = cvglmnetPredict(output.modelFits(iNeuron,1),X_full,params.lambdastring,'response');

    Yvars               = {};
    
    Yvars{1}            = Y_full;
    Yvars{2}            = Yh_full;
    
    for iC = 1:params.nShuffleCats
        idx             = ismember(output.x_label,params.shuffleVars{iC});
        X_temp          = X_full; 
        X_temp(:,~idx)  = 0;  
        Yh              = cvglmnetPredict(output.modelFits(iNeuron,1),X_temp,params.lambdastring,'response');
        Yvars{2+iC}     = Yh;
    end
%     Yvars{end+1}        = Y_full-Yh_full;
%     params.Yvarlabels{end+1} = 'Error';
    
    MOL_plotYYhat(temptrialData,params,Yvars,spikeData.cell_ID{iNeuron})
    
end

%% Compute variance explained over all individual trials:
% idx_time            = params.xtime>0 & params.xtime<1e6;
idx_time            = params.xtime>0 & params.xtime<=0.2e6;

cvR2_full           = NaN(nNeurons,1);
cvR2_full_time      = NaN(nNeurons,params.nTimebins);

cvR2_cat            = NaN(nNeurons,params.nShuffleCats);
cvR2_cat_time       = NaN(nNeurons,params.nShuffleCats,params.nTimebins);

cvR2_cat_un         = NaN(nNeurons,params.nShuffleCats);
cvR2_cat_un_time    = NaN(nNeurons,params.nShuffleCats,params.nTimebins);

fprintf('Computing cvR2 for neuron        \n');

% output.Yh_full      = NaN(size(output.y));
% output.Yh_cat       = NaN(size(output.y,1),size(output.y,2),params.nShuffleCats);

lastsesid = [];

for iNeuron = 1:nNeurons %loop over neurons
    fprintf(repmat('\b', 1, numel([num2str(iNeuron-1) num2str(nNeurons)])+2));
    fprintf('%d/%d\n',iNeuron,nNeurons);
    if ~isempty(output.modelFits(iNeuron,1).lambda)
        nTotalSpikeBins     = sum(~isnan(output.y(iNeuron,:))); %compute how many total spike bins there are from output struct
        nTrials             = nTotalSpikeBins / params.nTimebins; %how many trials iin this session
        %         idx_ses             = strcmp(sessionData.session_ID,spikeData.session_ID(iNeuron)); %get which session this is out of all ses predictors
        %         temptrialData       = MOL_getTempPerSes(sessionData.session_ID(idx_ses),trialData); %filter trialdata for this session
        temptrialData       = MOL_getTempPerSes(spikeData.session_ID(iNeuron),trialData);
        
        if ~(strcmp(spikeData.session_ID(iNeuron),lastsesid))
            load(fullfile(folderpath,sprintf('X_Full_GLM_%s.mat',spikeData.session_ID{iNeuron})))
            %  X_full              = squeeze(x(1,1:nTotalSpikeBins,1:params.nPredictors));
            lastsesid           = spikeData.session_ID(iNeuron);
        end
        
        %         X_full              = squeeze(output.x(idx_ses,1:nTotalSpikeBins,1:params.nPredictors)); %get predictor matrix for this session
        Y_full              = squeeze(output.y(iNeuron,1:nTotalSpikeBins))'; %get spike rates
        Yh_full             = cvglmnetPredict(output.modelFits(iNeuron,1),X_full,params.lambdastring,'response'); %cv predicted firing rates for full model
        idx_time_all        = repmat(idx_time,1,nTrials); %make an index of all time bins that are included for computation of cvR2
        cvR2_full(iNeuron,1) = 1 - var(Y_full(idx_time_all) - Yh_full(idx_time_all)) / var(Y_full(idx_time_all)); %compute R2

        %reshape and compute cvR2 for certain timebins:
        Y_r                 = reshape(Y_full,params.nTimebins,nTrials); 
        Yh_r                = reshape(Yh_full,params.nTimebins,nTrials);        
        cvR2_full_time(iNeuron,:) = 1 - var(Y_r - Yh_r,[],2) ./ var(Y_r,[],2); %compute R2
        
        for iC = 1:params.nShuffleCats
            %Predict firing rate when using only one predictor set, compute R2 and store value:
            idx                     = ismember(output.x_label,params.shuffleVars{iC});
            X_temp                  = X_full;
            X_temp(:,~idx)          = 0;
            Yh_cat                  = cvglmnetPredict(output.modelFits(iNeuron,1),X_temp,params.lambdastring,'response');
            cvR2_cat(iNeuron,iC)    = 1 - var(Y_full(idx_time_all) - Yh_cat(idx_time_all)) / var(Y_full(idx_time_all));
            
            Yh_cat_r                = reshape(Yh_cat,params.nTimebins,nTrials);
            cvR2_cat_time(iNeuron,iC,:) = 1 - var(Y_r - Yh_cat_r,[],2) ./ var(Y_r,[],2);
            
            %compute firing rate with circularly shuffled predictor and compute uniquely contributed R2:
            Yh_cat                  = cvglmnetPredict(output.shuffleFits(iNeuron,iC),X_full,params.lambdastring,'response');
            temp                    = 1 - var(Y_full(idx_time_all) - Yh_cat(idx_time_all)) / var(Y_full(idx_time_all));
            cvR2_cat_un(iNeuron,iC) = cvR2_full(iNeuron,1) - temp;
            
            Yh_cat_r                = reshape(Yh_cat,params.nTimebins,nTrials);
            temp                    = 1 - var(Y_r - Yh_cat_r,[],2) ./ var(Y_r,[],2);
            cvR2_cat_un_time(iNeuron,iC,:) = cvR2_full_time(iNeuron,:) - temp';
        end
    end
end

%% Correct for negative R2: 
% cvR2_full(cvR2_full<0)                  = 0;
% cvR2_full_time(cvR2_full_time<0)        = 0;
% cvR2_cat(cvR2_cat<0)                    = 0;
% cvR2_cat_time(cvR2_cat_time<0)          = 0;
% cvR2_cat_un(cvR2_cat_un<0)              = 0;
% cvR2_cat_un_time(cvR2_cat_un_time<0)    = 0;

%%
idx     = strcmp(spikeData.area,'V1');
fprintf('Variance explained across all V1 neurons (%2.3f, IQR %2.3f-%2.3f) (on trials excluding conflict trials) \n',...
    nanmedian(cvR2_full(idx,1)),prctile(cvR2_full(idx,1),[25 75]));

%% %remove animal 1012 with recordings in LM according to histology
spikeData.area(strcmp(spikeData.session_ID,'10122019041022') & strcmp(spikeData.area,'A1')) = deal({'LM'});
spikeData.area(strcmp(spikeData.session_ID,'10122019041126') & strcmp(spikeData.area,'A1')) = deal({'LM'});

%% Show cvR2 for each predictor subset and area:
params.nExperiments     = length(params.Experiments);
params.areas            = {'V1' 'A1'};% 'PPC' 'CG1'};
params.nAreas           = length(params.areas); 

figure; set(gcf,'units','normalized','Position',[0.2 0.43 0.22 0.3],'color','w'); hold all;

params.colors_splits = {[0.4 0.4 0.4] [0 0 0.8] [0.8 0 0] [0.1 0.7 0.3]};
handles = [];

for iArea = 1:params.nAreas
    idx     = strcmp(spikeData.area,params.areas{iArea});
    
%     idx                     = idx & ~(strcmp(spikeData.session_ID,'10122019041022') & strcmp(spikeData.area,'A1')); %remove animal 1012 with recordings in LM according to histology
%     idx                     = idx & ~(strcmp(spikeData.session_ID,'10122019041126') & strcmp(spikeData.area,'A1'));

    for iC = 1:params.nShuffleCats
        tmp     = cvR2_cat(idx,iC);
%         tmp     = cvR2_cat_un(idx,iC);
        
        handles(iC) = bar(iArea + iC/(params.nShuffleCats+1),nanmean(tmp),0.18,'k');
        set(handles(iC),'FaceColor',params.colors_splits{iC})
        errorbar(iArea + iC/(params.nShuffleCats+1),nanmean(tmp),nanstd(tmp)/sqrt(sum(idx)),'k','LineWidth',1,'CapSize',6);
    end
end

%statistics:
% Y               = cvR2_cat_un(:);
Y               = cvR2_cat(:);
X_var           = repmat(params.Yvarlabels(3:6),nNeurons,1);
X_var           = X_var(:);
X_area          = repmat(spikeData.area,1,4);
X_area          = X_area(:);
G_mou2          = repmat(G_mou,1,4);
G_mou2          = G_mou2(:);

idx_V1          = strcmp(X_area,'V1');
idx_AC          = strcmp(X_area,'A1');

tbl             = table(Y(idx_V1),X_var(idx_V1),G_mou2(idx_V1),'VariableNames',{'Y','Var','Mouse'}); %Create table for mixed model
lme             = fitlme(tbl,'Y~Var+(1|Mouse)'); %construct linear mixed effects model with fixed effect of temporal window (NO random intercept for mice as cohorts are tested)
stats           = dataset2table(anova(lme,'DFMethod','Satterthwaite')); %Perform ANOVA on model and output as matrix
fprintf('F(%d,%2.0f)=%1.2f, p=%1.2e; \n',stats{2,3},stats{2,4},stats{2,2},stats{2,5})

% 1.2000    1.4000    1.6000    1.8000
% 2.2000    2.4000    2.6000    2.8000

contrasts       = {[0 1 0 0] [0 0 1 0] [0 0 0 1] [0 1 -1 0] [0 1 0 -1] [0 0 1 -1]};
positions       = {[1.2 1.4] [1.2 1.6] [1.2 1.8] [1.4 1.6] [1.4 1.8] [1.6 1.8]};
contrastlabels  = {'T vs. V' 'T vs. A' 'T vs. M' 'V vs. A' 'V vs. M' 'A vs. M'};

fprintf('Posthoc comparison:\n')
for iC = 1:length(contrasts)
    [p,F,DF1,DF2] = coefTest(lme,contrasts{iC});
    fprintf('%s: F(%d,%d)=%1.1f, p=%1.2e; ',contrastlabels{iC},DF1,DF2,F,p)
    sigstar(positions{iC},p);
end

tbl             = table(Y(idx_AC),X_var(idx_AC),G_mou2(idx_AC),'VariableNames',{'Y','Var','Mouse'}); %Create table for mixed model
lme             = fitlme(tbl,'Y~Var+(1|Mouse)'); %construct linear mixed effects model with fixed effect of temporal window (NO random intercept for mice as cohorts are tested)
stats           = dataset2table(anova(lme,'DFMethod','Satterthwaite')); %Perform ANOVA on model and output as matrix
fprintf('F(%d,%2.0f)=%1.2f, p=%1.2e; \n',stats{2,3},stats{2,4},stats{2,2},stats{2,5})

fprintf('Posthoc comparison:\n')
for iC = 1:length(contrasts)
    [p,F,DF1,DF2] = coefTest(lme,contrasts{iC});
    fprintf('%s: F(%d,%d)=%1.1f, p=%1.2e; ',contrastlabels{iC},DF1,DF2,F,p)
    sigstar(positions{iC}+1,p);
end

for iC = 1:4
    idx = strcmp(X_var,params.Yvarlabels(2+iC));
    tbl             = table(Y(idx),X_area(idx),G_mou2(idx),'VariableNames',{'Y','Area','Mouse'}); %Create table for mixed model
    lme             = fitlme(tbl,'Y~Area+(1|Mouse)'); %construct linear mixed effects model with fixed effect of temporal window (NO random intercept for mice as cohorts are tested)
    stats           = dataset2table(anova(lme,'DFMethod','Satterthwaite')); %Perform ANOVA on model and output as matrix
    fprintf('\n%s: ',params.Yvarlabels{2+iC})
    fprintf('F(%d,%2.0f)=%1.1f, p=%1.2e;',stats{2,3},stats{2,4},stats{2,2},stats{2,5})
    sigstar([1 2] + iC/(params.nShuffleCats+1),stats{2,5});
end

set(gca,'XTick',1.5:1:5.5,'XTickLabels',params.areas,'YTick',0:0.01:0.2)
xlim([1 1+params.nAreas])
% ylim([0 0.06])
legend(handles(1:4),{'Trial' 'Visual' 'Audio' 'Motor'},'Location','NorthWest'); legend boxoff;

% export_fig(fullfile(params.savedir,sprintf('Bar_cvR2_V1_AC')),'-eps','-nocrop')

%% Show cvR2 for each predictor subset and area for each cohort separately:
params.nExperiments     = length(params.Experiments);
params.areas            = {'V1' 'A1'};
params.nAreas           = length(params.areas); 

figure; set(gcf,'units','normalized','Position',[0.2 0.43 0.45 0.3],'color','w'); hold all;

params.colors_splits = {[0.4 0.4 0.4] [0 0 0.8] [0.8 0 0] [0.1 0.7 0.3]};
handles = [];

for iExp = 1:3
    subplot(1,3,iExp); hold all;
    idx_exp                 = ismember(spikeData.session_ID,sessionData.session_ID(strcmp(sessionData.Experiment,params.Experiments(iExp))));
    
    for iArea = 1:params.nAreas
        idx_area = strcmp(spikeData.area,params.areas{iArea});
        idx     = idx_exp & idx_area;
        for iC = 1:params.nShuffleCats
%             tmp     = cvR2_cat(idx,iC);
            tmp     = cvR2_cat_un(idx,iC);
            handles(iC) = bar(iArea + iC/(params.nShuffleCats+1),nanmean(tmp),0.18,'k');
            set(handles(iC),'FaceColor',params.colors_splits{iC})
            errorbar(iArea + iC/(params.nShuffleCats+1),nanmean(tmp),nanstd(tmp)/sqrt(sum(idx)),'k','LineWidth',1,'CapSize',6);
        end
    end
    set(gca,'XTick',1.5:1:5.5,'XTickLabels',params.areas,'YTick',0:0.01:0.2)
    xlim([1 1+params.nAreas])
    ylim([0 0.08])
    title(params.ExperimentLabels(iExp),'FontSize',10)
%     legend(handles(1:4),{'Trial' 'Visual' 'Audio' 'Motor'},'Location','NorthWest'); legend boxoff;
    if iExp == 1
        legend(handles(1:4),{'Trial' 'Visual' 'Audio' 'Motor'},'Location','best'); legend boxoff;
    end
end
% export_fig(fullfile(params.savedir,sprintf('Bar_cvR2_V1_AC_3Cohorts')),'-eps','-nocrop')

%% Show cvR2 for each predictor subset and area for each cohort separately:
params.areas            = {'V1' 'A1'};

meantoplot = NaN(3,2);
for iExp = 1:3
    idx_exp                 = ismember(spikeData.session_ID,sessionData.session_ID(strcmp(sessionData.Experiment,params.Experiments(iExp))));
    
    idx_area                = strcmp(spikeData.area,'V1');
    idx                     = idx_exp & idx_area;
    meantoplot(iExp,1)      = nanmean(cvR2_cat(idx,3));
    meantoplot(iExp,2)      = nanmean(cvR2_cat(idx,4));
end

figure; set(gcf,'units','normalized','Position',[0.2 0.43 0.08 0.15],'color','w'); hold all;

h = bar(meantoplot,0.7,'stacked'); %plot bars as stacked:
for iExp = 1:3 %add text of percentage explained:
    for iVar = 1:2
        text(iExp-0.3,h(iVar).YData(iExp),sprintf('%2.1f%%',(meantoplot(iExp,iVar)/sum(meantoplot(iExp,:)))*100),'FontSize',6)
    end
end
h(1).FaceColor = params.colors_splits{3};
h(2).FaceColor = params.colors_splits{4};
xlim([0.5 3.5])
set(gca,'XTick',1:3,'XTickLabels',params.ExperimentLabels)

export_fig(fullfile(params.savedir,sprintf('Bar_cvR2_AuResp_Perc_3Cohorts')),'-eps','-nocrop')

%% Show explained variance for each category of predictors over time:
% !DEPRECATED
params.labels_cats      = {'Full' 'Trial' 'Vis' 'Aud' 'Motor'};
params.colors_cats      = {[0 0 0] [0.5 0.2 0.3] [0 0 0.8] [0.8 0 0] [0.1 0.7 0.3]};

% iExp = 2;
% idx_exp                 = ismember(spikeData.session_ID,sessionData.session_ID(strcmp(sessionData.Experiment,params.Experiments(iExp))));
idx_exp                 = ismember(spikeData.session_ID,sessionData.session_ID(ismember(sessionData.Experiment,params.Experiments)));

idx_area                = strcmp(spikeData.area,'V1');
idx                     = idx_exp & idx_area;

figure; set(gcf,'units','normalized','Position',[0.05 0.59 0.37 0.25],'color','w'); hold all;
subplot(1,2,1); hold all;
meantoplot      = nanmean(cvR2_full_time(idx,:),1);
errortoplot     = nanstd(cvR2_full_time(idx,:),[],1) / sqrt(sum(idx));
plot(params.xtime,meantoplot,'Color',params.colors_cats{1},'LineWidth',0.5);

for iC = 1:4
    meantoplot      = squeeze(nanmean(cvR2_cat_time(idx,iC,:),1));
%     meantoplot      = squeeze(nanmean(cvR2_cat_un_time(idx,iC,:),1));
    errortoplot     = squeeze(nanstd(cvR2_cat_time(idx,iC,:),[],1) / sqrt(sum(idx)));
    plot(params.xtime,meantoplot,'Color',params.colors_cats{iC+1},'LineWidth',0.5)
end

% meantoplot      = squeeze(nanmean(nansum(cvR2_cat_time(idx,:,:),2),1));
% plot(params.xtime,meantoplot,'Color',[1 0 1],'LineWidth',0.5)

xlim([params.t_pre params.t_post])
ylim([0 0.1])
plot([0 0],get(gca,'ylim'),'k:','LineWidth',0.5)
set(gca,'XTick',-2e6:0.2e6:2e6,'XTickLabels',(-2e6:0.2e6:2e6)*1e-3,'YTick',get(gca,'ylim'));
xlabel('Time (ms)')
title('V1')

idx_area = strcmp(spikeData.area,'A1');
idx         = idx_exp & idx_area;

subplot(1,2,2); hold all;
meantoplot      = nanmean(cvR2_full_time(idx,:),1);
errortoplot     = nanstd(cvR2_full_time(idx,:),[],1) / sqrt(sum(idx));
plot(params.xtime,meantoplot,'Color',params.colors_cats{1},'LineWidth',0.5);

for iC = 1:4
    meantoplot      = squeeze(nanmean(cvR2_cat_time(idx,iC,:),1));
%     meantoplot      = squeeze(nanmean(cvR2_cat_un_time(idx,iC,:),1));
    errortoplot     = squeeze(nanstd(cvR2_cat_time(idx,iC,:),[],1) / sqrt(sum(idx)));
    plot(params.xtime,meantoplot,'Color',params.colors_cats{iC+1},'LineWidth',0.5)
end
xlim([params.t_pre params.t_post])
ylim([0 0.1])
plot([0 0],get(gca,'ylim'),'k:','LineWidth',0.5)
set(gca,'XTick',-2e6:0.2e6:2e6,'XTickLabels',(-2e6:0.2e6:2e6)*1e-3,'YTick',get(gca,'ylim'));
xlabel('Time (ms)')
title('AC')
legend(params.labels_cats,'FontSize',12); legend boxoff;

export_fig(fullfile(params.savedir,sprintf('cvR2_overTime_V1_AC')),'-eps','-nocrop')

%% Show explained variance across neurons as bars:
params.labels_cats      = {'Raw' 'Full' 'Trial' 'Vis' 'Aud' 'Motor'};
params.colors_cats      = {[0 0 0] [0.5 0.5 0.5] [0.5 0.2 0.3] [0 0 0.8] [0.8 0 0] [0.1 0.7 0.3]};

idx_exp                 = true(size(spikeData.session_ID)); %all cohorts

idx_area                = strcmp(spikeData.area,'V1');
idx                     = idx_exp & idx_area;

%Sort by aud/vis/motor
% temp                    = cvR2_cat(idx,:);
temp                    = cvR2_cat_un(idx,:);
sortval                 = (-(temp(:,2) / max(temp(:,2))) + (temp(:,4) / max(temp(:,4)))) ./ (temp(:,3) / max(temp(:,3)));
[~,sortidx]             = sort(sortval);
temp                    = temp(sortidx,:); %actual sorting based on derived index

figure; set(gcf,'units','normalized','Position',[0.05 0.59 0.16 0.2],'color','w'); hold all;
for iC = 2:4 %for each of the categories (except trialData)
    subplot(3,1,iC-1); hold all;
    h = bar(1:sum(idx),temp(:,iC),1.5,'k');
    ylim(round(get(gca,'ylim'),1));
    set(gca,'XTick',[],'YTick',get(gca,'ylim'));
    ylabel('cvR2')
%     title(params.labels_cats(iC+2))
    text(sum(idx)/2.5,0.06,params.labels_cats(iC+2),'FontSize',8)
end

export_fig(fullfile(params.savedir,sprintf('Bars_cvR2_AV_3cohorts_confexcl')),'-eps','-nocrop')

comparisons = [2 4; 2 3; 3 4];
for i = 1:3
    %statistics:
%     tbl             = table(cvR2_cat(idx,comparisons(i,1)),cvR2_cat(idx,comparisons(i,2)),G_mou(idx),'VariableNames',{'Y','X','Mouse'}); %Create table for mixed model
    tbl             = table(cvR2_cat_un(idx,comparisons(i,1)),cvR2_cat_un(idx,comparisons(i,2)),G_mou(idx),'VariableNames',{'Y','X','Mouse'}); %Create table for mixed model
    lme             = fitlme(tbl,'Y~X+(1|Mouse)'); %construct linear mixed effects model with fixed effect of temporal window (NO random intercept for mice as cohorts are tested)
    stats           = dataset2table(anova(lme,'DFMethod','Satterthwaite')); %Perform ANOVA on model and output as matrix
    r               = sqrt(lme.Rsquared.Ordinary);
    
    fprintf('Correlation between %s and %s \n',params.labels_cats{comparisons(i,1)+2},params.labels_cats{comparisons(i,2)+2})
    if stats{2,5}<0.001
        fprintf('r=%1.2f, F(%d,%2.0f)=%1.2f, p=%1.2e; \n',r,stats{2,3},stats{2,4},stats{2,2},stats{2,5})
    else
        fprintf('r=%1.2f, F(%d,%2.0f)=%1.2f, p=%1.3f; \n',r,stats{2,3},stats{2,4},stats{2,2},stats{2,5})
    end
end

%% Show explained variance across neurons as scatter based on rank: (stringer et al. 2019 S13)
idx_exp                 = true(size(spikeData.session_ID)); %all cohorts
idx_area                = strcmp(spikeData.area,'V1');
idx                     = idx_exp & idx_area;

% temp = cvR2_cat_un(:,:);
% temp = temp + randn(size(temp))*1e-6;

I = NaN(sum(idx),4);
for iC = 1:4 %for each of the categories find the ranks
    I(:,iC) = tiedrank(cvR2_cat(idx,iC));
%     I(:,iC) = tiedrank(temp(idx,iC));
end

figure; set(gcf,'units','normalized','Position',[0.05 0.6 0.37 0.17],'color','w'); hold all;
comparisons = [2 4; 2 3; 3 4];
for iC = 1:3 %for each of the category combinations (except trialData)
    subplot(1,3,iC); hold all;
    h = scatter(I(:,comparisons(iC,1)),I(:,comparisons(iC,2)),10,'k.');
    ylim([0 sum(idx)]);
    xlim([0 sum(idx)]);
    set(gca,'XTick',get(gca,'xlim'),'YTick',get(gca,'ylim'));
    xlabel(sprintf('from %s (rank)',params.labels_cats{comparisons(iC,1)+2}))
    ylabel(sprintf('from %s (rank)',params.labels_cats{comparisons(iC,2)+2}))
end

export_fig(fullfile(params.savedir,sprintf('Scatter_cvR2_AV_StringerS13')),'-eps','-nocrop')

comparisons = [2 4; 2 3; 3 4];
for i = 1:3
    %statistics:
    [r,p] = corr(cvR2_cat(idx,comparisons(i,1)), cvR2_cat(idx,comparisons(i,2)), 'type', 'Spearman','rows' ,'complete');
%     [r,p] = corr(temp(idx,comparisons(i,1)), temp(idx,comparisons(i,2)), 'type', 'Spearman','rows' ,'complete');
    fprintf('Spearman rank correlation between %s and %s \n',params.labels_cats{comparisons(i,1)+2},params.labels_cats{comparisons(i,2)+2})
    if p<0.001
        fprintf('r=%1.2f, p=%1.2e; \n',r,p)
    else
        fprintf('r=%1.2f, p=%1.3f; \n',r,p)
    end
end

%% Computing significance of EV versus shuffling firing rate and prediction:
nShuffles       = 1000;
alphathr        = 0.05;
nNeurons        = length(spikeData.session_ID);

idx_time        = params.xtime>0 & params.xtime<0.2e6;

% temp_var_expl_splits = var_expl_splits(expidx,:);
cvR2_shuf         = NaN(nNeurons,nShuffles);
cvR2_cat_shuf     = NaN(nNeurons,params.nShuffleCats,nShuffles);

for iNeuron = 1:nNeurons %loop over neurons
    fprintf(repmat('\b', 1, numel([num2str(iNeuron-1) num2str(nNeurons)])+2));
    fprintf('%d/%d\n',iNeuron,nNeurons);
    if ~isempty(output.modelFits(iNeuron,1).lambda)
        nTotalSpikeBins     = sum(~isnan(output.y(iNeuron,:))); %compute how many total spike bins there are from output struct
        nTrials             = nTotalSpikeBins / params.nTimebins; %how many trials iin this session
        temptrialData       = MOL_getTempPerSes(spikeData.session_ID(iNeuron),trialData);
        
        if ~(strcmp(spikeData.session_ID(iNeuron),lastsesid))
            load(fullfile(folderpath,sprintf('X_Full_GLM_%s.mat',spikeData.session_ID{iNeuron})))
            lastsesid           = spikeData.session_ID(iNeuron);
        end
        
        idx_time_all        = repmat(idx_time,1,nTrials); %make an index of all time bins that are included for computation of cvR2

%         Yh_full             = cvglmnetPredict(output.modelFits(iNeuron,1),X_full,params.lambdastring,'response'); %cv predicted firing rates for full model
% 
%         for iS = 1:nShuffles
%             Y_full_shuf                = squeeze(output.y(iNeuron,1:nTotalSpikeBins))'; %get spike rates
%             Y_full_shuf                = reshape(Y_full_shuf,params.nTimebins,nTrials); %reshape to time by trial
%             Y_full_shuf                = Y_full_shuf(:,randperm(nTrials)); %permute these trials, keep NaNs for higher trial numbers the same
%             Y_full_shuf                = reshape(Y_full_shuf,nTotalSpikeBins,1); %reshape again to vector
%             
%             cvR2_shuf(iNeuron,iS) = 1 - var(Y_full_shuf(idx_time_all) - Yh_full(idx_time_all)) / var(Y_full_shuf(idx_time_all)); %compute R2
%         end
        
        for iC = 1:params.nShuffleCats
            %Predict firing rate when using only one predictor set, compute R2 and store value:
            idx                     = ismember(output.x_label,params.shuffleVars{iC});
            X_temp                  = X_full;
            X_temp(:,~idx)          = 0;
            Yh_cat                  = cvglmnetPredict(output.modelFits(iNeuron,1),X_temp,params.lambdastring,'response');
            
            for iS = 1:nShuffles
                Y_full_shuf                = squeeze(output.y(iNeuron,1:nTotalSpikeBins))'; %get spike rates
                Y_full_shuf                = reshape(Y_full_shuf,params.nTimebins,nTrials); %reshape to time by trial
                Y_full_shuf                = Y_full_shuf(:,randperm(nTrials)); %permute these trials, keep NaNs for higher trial numbers the same
                Y_full_shuf                = reshape(Y_full_shuf,nTotalSpikeBins,1); %reshape again to vector
                
                cvR2_cat_shuf(iNeuron,iC,iS) = 1 - var(Y_full_shuf(idx_time_all) - Yh_cat(idx_time_all)) / var(Y_full_shuf(idx_time_all)); %compute R2
            end
        end
    end
end

%% Determine significance:
cvR2_cat_sig     = cvR2_cat > prctile(cvR2_cat_shuf,99.9,3);

%% Reconstruct index to get the appropriate data:
idx_exp                 = true(size(spikeData.session_ID)); %all cohorts
idx_area                = strcmp(spikeData.area,'V1');
idx                     = idx_exp & idx_area;

%% 
nShuffles       = 1000;

pairs           = {2 3; 2 4; 3 4;};
nNeurons        = sum(idx);

mat_overlap_shuf     = NaN(3,nShuffles);

for iS = 1:nShuffles
    temp                    = cvR2_cat_sig(idx,:);
    for iC = 1:4
        temp(:,iC)           = temp(randperm(nNeurons),iC);
    end
    
    for iP = 1:3
        mat_overlap_shuf(iP,iS)     = sum(temp(:,pairs{iP,1}) & temp(:,pairs{iP,2})) / sum(idx);
    end
    
end

for iP = 1:3
    mat_overlap_orig(iP)     = sum(cvR2_cat_sig(idx,pairs{iP,1}) & cvR2_cat_sig(idx,pairs{iP,2})) / sum(idx);
end

for iP = 1:3
    p = sum(mat_overlap_orig(iP)>mat_overlap_shuf(iP,:))/1000;
    fprintf('Overlap between %s and %s: %1.3f \n',params.Yvarlabels{pairs{iP,1}+2},params.Yvarlabels{pairs{iP,2}+2},p)
end

%% Compute predicted firing rate response for trial types based on different subpredictors:

fprintf('Computing firing rate response based on subpredictors for neuron        \n');

params.nSplits      = 4;
ratemat             = NaN(nNeurons,params.nShuffleCats+2,params.nSplits,params.nTimebins);

%Store for each split whether the 0-200 ms response is different from baseline. Third dimension different predictions:
%Types are: original data, full model, model without movement, original data - movement prediction
signmat             = false(nNeurons,params.nSplits,4); 

params.minTrialCond = 10;

params.twin_baseline_start = -500e3;
params.twin_baseline_stop = 0;
params.twin_resp_start = 0;
params.twin_resp_stop = 200e3;

params.ttestalpha = 0.025;

lastsesid   = [];
params.lambdastring = 'lambda_min';

for iNeuron = 1:nNeurons %loop over neurons
    fprintf(repmat('\b', 1, numel([num2str(iNeuron-1) num2str(nNeurons)])+2));
    fprintf('%d/%d\n',iNeuron,nNeurons);
    if ~isempty(output.modelFits(iNeuron,1).lambda)
        
        nTotalSpikeBins     = sum(~isnan(output.y(iNeuron,:)));
        nTrials             = nTotalSpikeBins / params.nTimebins;
        
        temptrialData       = MOL_getTempPerSes(spikeData.session_ID(iNeuron),trialData);

        splits              = {};
        splits{1}           = ismember(temptrialData.trialType,{'X'}) & temptrialData.visualOriChangeNorm==2;
        splits{2}           = ismember(temptrialData.trialType,{'X'}) & temptrialData.visualOriChangeNorm==3;
        
        splits{3}           = ismember(temptrialData.trialType,{'Y'}) & temptrialData.audioFreqChangeNorm==2;
        splits{4}           = ismember(temptrialData.trialType,{'Y'}) & temptrialData.audioFreqChangeNorm==3;
        
        if ~(strcmp(spikeData.session_ID(iNeuron),lastsesid))
            load(fullfile(folderpath,sprintf('X_Full_GLM_%s.mat',spikeData.session_ID{iNeuron})))
%             X_full              = squeeze(x(1,1:nTotalSpikeBins,1:params.nPredictors));
            lastsesid           = spikeData.session_ID(iNeuron);
        end
        
        %         X_full              = squeeze(output.x(idx_ses,1:nTotalSpikeBins,1:params.nPredictors));
        Y_full              = squeeze(output.y(iNeuron,1:nTotalSpikeBins))';
        
        Y_r                 = reshape(Y_full,params.nTimebins,nTrials);
        
        for iSplit = 1:params.nSplits
            ratemat(iNeuron,1,iSplit,:)     = nanmean(Y_r(:,splits{iSplit}),2);
        end
        
        Yh_full             = cvglmnetPredict(output.modelFits(iNeuron,1),X_full,params.lambdastring,'response');
        Yh_r                = reshape(Yh_full,params.nTimebins,nTrials);
        
        for iSplit = 1:params.nSplits
            ratemat(iNeuron,2,iSplit,:)     = nanmean(Yh_r(:,splits{iSplit}),2);
        end
        
        for iC = 1:params.nShuffleCats
            %Predict firing rate when using only one predictor set, compute R2 and store value:
            idx                     = ismember(output.x_label,params.shuffleVars{iC});
            X_temp                  = X_full;
            X_temp(:,~idx)          = 0;
            Yh_cat                  = cvglmnetPredict(output.modelFits(iNeuron,1),X_temp,params.lambdastring,'response');
            Yh_cat_r                = reshape(Yh_cat,params.nTimebins,nTrials);

            for iSplit = 1:params.nSplits
                ratemat(iNeuron,iC+2,iSplit,:)     = nanmean(Yh_cat_r(:,splits{iSplit}),2);
            end
        end
        
        %Compute fracion of significantly responsive neurons on original as well as when regressing out movement variability:
        %First predict firing rate when using full model expect for motor-predictors:
        iC = 4;
        idx                     = ismember(output.x_label,params.shuffleVars{iC});
        X_temp                  = X_full;
        X_temp(:,idx)           = 0;
        Yh_nomot                = cvglmnetPredict(output.modelFits(iNeuron,1),X_temp,params.lambdastring,'response');
        Yh_nomot_r              = reshape(Yh_nomot,params.nTimebins,nTrials);
        
        X_temp                  = X_full;
        X_temp(:,~idx)          = 0;
        Yh_mot                  = cvglmnetPredict(output.modelFits(iNeuron,1),X_temp,params.lambdastring,'response');
        Yh_mot_r                = reshape(Yh_mot,params.nTimebins,nTrials);

        Yh_submot_r             = Y_r - Yh_mot_r;
        
        %Compute baseline activity using original or predicted data without movement:
        varstrings = {'Y_r' 'Yh_r' 'Yh_nomot_r' 'Yh_submot_r'};
        
        for iVar = 1:size(signmat,3)
            for iSplit = 1:params.nSplits %Store the mean response for each of these splits
                if sum(splits{iSplit})>=params.minTrialCond
                    eval(sprintf('tempdata = %s'';',varstrings{iVar}))
                    bsl                                     = nanmean(tempdata(splits{iSplit},params.xtime>params.twin_baseline_start & params.xtime<params.twin_baseline_stop),2);
                    resp                                    = nanmean(tempdata(splits{iSplit},params.xtime>params.twin_resp_start & params.xtime<params.twin_resp_stop),2);
                    [~,signmat(iNeuron,iSplit,iVar)]        = signrank(bsl,resp,'alpha',params.ttestalpha);
                end
            end
        end
    end
end

%% Show mean firing rate response to trial types with predictions based on different subsets. 

params.labels_splits    = {'Vthr' 'Vmax' 'Athr' 'Amax'};
params.labels_cats      = {'Raw' 'Full' 'Trial' 'Vis' 'Aud' 'Motor'};
params.colors_cats      = {[0 0 0] [0.5 0.5 0.5] [0.5 0.2 0.3] [0 0 0.8] [0.8 0 0] [0.1 0.7 0.3]};

idx_exp                 = true(size(spikeData.session_ID));

idx_area                = strcmp(spikeData.area,'V1');
idx                     = idx_exp & idx_area;

figure; set(gcf,'units','normalized','Position',[0.25 0.5 0.45 0.28],'color','w'); hold all;

for iSplit = 1:params.nSplits
    subplot(2,params.nSplits,iSplit); hold all;
    trialoffset = squeeze(nanmean(nanmean(diff(ratemat(idx,[1 3],iSplit,params.xtime<0),[],2),1),4));
    
    for iC = [1 2 4 5 6]
        tmp             = squeeze(nanmean(ratemat(idx,iC,iSplit,:),1));
        if iC>2 %correct for trial prediction in display:
            tmp = tmp + trialoffset;
        end
        plot(params.xtime,tmp,'Color',params.colors_cats{iC},'LineWidth',0.5)
    end
    xlim([params.t_pre params.t_post])
    ylim([4.5 9])
    plot([0 0],get(gca,'ylim'),'k:','LineWidth',0.5)
    if iSplit ==4
        set(gca,'XTick',-2e6:0.2e6:2e6,'XTickLabels',(-2e6:0.2e6:2e6)*1e-3, 'YTick',get(gca,'ylim'));
    else
        set(gca,'XTick',-2e6:0.2e6:2e6,'XTickLabels',[],'YTick',[]);
    end
    if iSplit==params.nSplits
        legend(params.labels_cats([1 2 4 5 6]),'FontSize',8);
        legend boxoff
    end
    title(params.labels_splits{iSplit},'FontSize',9)
end

idx_area    = strcmp(spikeData.area,'A1');
idx         = idx_exp & idx_area;

% figure; set(gcf,'units','normalized','Position',[0.05 0.19 0.27 0.28],'color','w'); hold all;
for iSplit = 1:params.nSplits
%     subplot(2,params.nSplits/2,iSplit); hold all;
    subplot(2,params.nSplits,iSplit+params.nSplits); hold all;
    
    trialoffset = squeeze(nanmean(nanmean(diff(ratemat(idx,[1 3],iSplit,params.xtime<0),[],2),1),4));
    
    for iC = [1 2 4 5 6]
        tmp             = squeeze(nanmean(ratemat(idx,iC,iSplit,:),1));
        if iC>2 %correct for trial prediction in display:
            tmp = tmp - trialoffset;
        end
        plot(params.xtime,tmp,'Color',params.colors_cats{iC},'LineWidth',0.5)
    end
    
    xlim([params.t_pre params.t_post])
    ylim([3.5 7])
    plot([0 0],get(gca,'ylim'),'k:','LineWidth',0.5)
    if iSplit == 4
        set(gca,'XTick',-2e6:0.2e6:2e6,'XTickLabels',(-2e6:0.2e6:2e6)*1e-3, 'YTick',get(gca,'ylim'));
    else
        set(gca,'XTick',-2e6:0.2e6:2e6,'XTickLabels',[],'YTick',[]);
    end
    
    title(params.labels_splits{iSplit},'FontSize',9)
end

export_fig(fullfile(params.savedir,sprintf('Rate_overTime_V1_AC_4trialtypes')),'-eps','-nocrop')

%% For the three cohorts separately: Show mean firing rate response to trial types with predictions based on different subsets. 

params.labels_splits    = {'Vthr' 'Vmax' 'Athr' 'Amax'};
params.labels_cats      = {'Raw' 'Full' 'Trial' 'Vis' 'Aud' 'Motor'};
params.colors_cats      = {[0 0 0] [0.5 0.5 0.5] [0.5 0.2 0.3] [0 0 0.8] [0.8 0 0] [0.1 0.7 0.3]};

idx_area                = strcmp(spikeData.area,'V1');

figure; set(gcf,'units','normalized','Position',[0.25 0.5 0.45 0.28],'color','w'); hold all;

for iExp = 1:3
    
    idx_exp                 = ismember(spikeData.session_ID,sessionData.session_ID(strcmp(sessionData.Experiment,params.Experiments(iExp))));
    % idx_exp                 = ismember(spikeData.session_ID,sessionData.session_ID(ismember(sessionData.Experiment,params.Experiments)));
    
    idx                     = idx_exp & idx_area;
    
%     for iSplit = 1:params.nSplits
    for iSplit = [3 4]
        subplot(3,params.nSplits,(iExp-1)*params.nSplits + iSplit); hold all;
%         trialoffset = squeeze(nanmean(nanmean(diff(ratemat(idx,[1 3],iSplit,params.xtime<0),[],2),1),4));
        trialoffset = squeeze(nanmean(nanmean(diff(ratemat(idx,[5 3],iSplit,params.xtime<0),[],2),1),4));
        for iC = [1 2 4 5 6]
            tmp             = squeeze(nanmean(ratemat(idx,iC,iSplit,:),1));
            if iC>3 %correct for trial prediction in display:
                tmp = tmp + trialoffset;
            end
            plot(params.xtime,tmp,'Color',params.colors_cats{iC},'LineWidth',0.5)
        end
        xlim([params.t_pre+0.1e6 params.t_post])
        temp = get(gca,'ylim');
        ylim([floor(temp(1)) ceil(temp(2))]);
%             ylim([3 9])
        plot([0 0],get(gca,'ylim'),'k:','LineWidth',0.5)
%         if iSplit == 3
            set(gca,'XTick',-2e6:0.2e6:2e6,'XTickLabels',(-2e6:0.2e6:2e6)*1e-3, 'YTick',get(gca,'ylim'));
%         else
%             set(gca,'XTick',-2e6:0.2e6:2e6,'XTickLabels',[],'YTick',get(gca,'ylim'),'YTickLabels',[]);
%         end
        if iSplit==params.nSplits && iExp==3
            legend(params.labels_cats([1 2 4 5 6]),'FontSize',8);
            legend boxoff
        end
        title(params.labels_splits{iSplit},'FontSize',9)
    end
end

export_fig(fullfile(params.savedir,sprintf('Rate_overTime_V1_3cohorts_4trialtypes')),'-eps','-nocrop')

%% For the three cohorts separately: show auditory response

idx_area                = strcmp(spikeData.area,'V1');

response                = squeeze(nanmean(ratemat(:,:,:,params.xtime>0 & params.xtime<=200e3),4)- nanmean(ratemat(:,:,:,params.xtime<0),4));

figure; set(gcf,'units','normalized','Position',[0.25 0.5 0.2 0.18],'color','w'); hold all;

handles = [];

subplot(1,2,1); hold all;
datatotest = squeeze(nanmean(response(idx_area,5,[3 4]),3));

for iExp = 1:3
    idx_exp                 = ismember(spikeData.session_ID,sessionData.session_ID(strcmp(sessionData.Experiment,params.Experiments(iExp))));
    idx                     = idx_area & idx_exp;
    datatoplot              = squeeze(nanmean(response(idx,5,[3 4]),3));
    h = bar(iExp,nanmean(datatoplot),'FaceColor',params.colors_experiments{iExp});
    errorbar(iExp,nanmean(datatoplot),nanstd(datatoplot)/sqrt(sum(idx_exp)),'k','LineWidth',1)
end
ylim([0 0.8])

tbl             = table(datatotest,X_coh(idx_area),G_mou(idx_area),'VariableNames',{'Y','Var','Mouse'}); %Create table for mixed model
lme             = fitlme(tbl,'Y~Var'); %construct linear mixed effects model with fixed effect of temporal window (NO random intercept for mice as cohorts are tested)
stats           = dataset2table(anova(lme,'DFMethod','Satterthwaite')); %Perform ANOVA on model and output as matrix
fprintf('F(%d,%2.0f)=%1.2f, p=%1.2e; \n',stats{2,3},stats{2,4},stats{2,2},stats{2,5})
sigstar([1 3],stats{2,5})
set(gca,'XTick',1:3,'XTickLabels',params.ExperimentLabels)

subplot(1,2,2); hold all;
datatotest = squeeze(nanmean(response(idx_area,6,[3 4]),3));

for iExp = 1:3
    idx_exp                 = ismember(spikeData.session_ID,sessionData.session_ID(strcmp(sessionData.Experiment,params.Experiments(iExp))));
    idx                     = idx_area & idx_exp;
    datatoplot              = squeeze(nanmean(response(idx,6,[3 4]),3));
    h = bar(iExp,nanmean(datatoplot),'FaceColor',params.colors_experiments{iExp});
    errorbar(iExp,nanmean(datatoplot),nanstd(datatoplot)/sqrt(sum(idx_exp)),'k','LineWidth',1)
end
ylim([0 0.8])
tbl             = table(datatotest,X_coh(idx_area),G_mou(idx_area),'VariableNames',{'Y','Var','Mouse'}); %Create table for mixed model
lme             = fitlme(tbl,'Y~Var'); %construct linear mixed effects model with fixed effect of temporal window (NO random intercept for mice as cohorts are tested)
stats           = dataset2table(anova(lme,'DFMethod','Satterthwaite')); %Perform ANOVA on model and output as matrix
fprintf('F(%d,%2.0f)=%1.2f, p=%1.2e; \n',stats{2,3},stats{2,4},stats{2,2},stats{2,5})

contrasts       = {[0 1 0] [0 0 1] [0 1 -1]};
positions       = {[1 3] [1 2] [2 3]};
contrastlabels  = {'NE vs. MST' 'NE vs. UST' 'UST vs. MST'};

fprintf('Posthoc comparison:\n')
for iC = 1:length(contrasts)
    [p,F,DF1,DF2] = coefTest(lme,contrasts{iC});
    fprintf('%s: F(%d,%d)=%1.1f, p=%1.2e; \n',contrastlabels{iC},DF1,DF2,F,p)
    sigstar(positions{iC},p);
end

ylabel('Evoked response (sp/s)')
set(gca,'XTick',1:3,'XTickLabels',params.ExperimentLabels)

export_fig(fullfile(params.savedir,sprintf('Bar_Rate_V1_3Cohorts')),'-eps','-nocrop')


%% Show explained variance for each category of predictors over time across neurons as a heatmap:

% iExp = 3;
% idx_exp                 = ismember(spikeData.session_ID,sessionData.session_ID(strcmp(sessionData.Experiment,params.Experiments(iExp))));
idx_exp                 = ismember(spikeData.session_ID,sessionData.session_ID(ismember(sessionData.Experiment,params.Experiments)));

idx_area                = strcmp(spikeData.area,'V1');
idx                     = idx_exp & idx_area;

%Zscore the firing rates:
temp                    = ratemat - repmat(nanmean(ratemat(:,:,:,params.xtime<0),4),1,1,1,params.nTimebins);
temp2                   = ratemat(:,:,:,params.xtime<0);
temp                    = temp ./ repmat(nanstd(temp2(:)),nNeurons,params.nShuffleCats+2,params.nSplits,params.nTimebins); %#ok<*NANSTD>

%Make heatmap for visual trials::
% idx_trialtypes          = [3 4];
idx_ttypes_vis          = [1 2];
idx_ttypes_aud          = [3 4];

heatmat_vis             = ones(sum(idx),params.nTimebins,3); %init
temp3                   = nanmean(temp(idx,4,idx_ttypes_vis,:),3);
heatmat_vis(:,:,1)      = heatmat_vis(:,:,1) - squeeze(nanmean(temp(idx,4,idx_ttypes_vis,:),3)) / prctile(temp3(:),98);
heatmat_vis(:,:,2)      = heatmat_vis(:,:,2) - squeeze(nanmean(temp(idx,4,idx_ttypes_vis,:),3)) / prctile(temp3(:),98);
temp4                   = nanmean(temp(idx,6,idx_ttypes_vis,:),3);
heatmat_vis(:,:,1)      = heatmat_vis(:,:,1) - squeeze(nanmean(temp(idx,6,idx_ttypes_vis,:),3)) / prctile(temp4(:),98);
heatmat_vis(:,:,3)      = heatmat_vis(:,:,3) - squeeze(nanmean(temp(idx,6,idx_ttypes_vis,:),3)) / prctile(temp4(:),98);

%sort:
% sortval(:,1)            = nanmean(nanmean(temp(idx,4,idx_trialtypes,params.xtime>0 & params.xtime<0.2e6),4),3);
% sortval(:,2)            = nanmean(nanmean(temp(idx,6,idx_trialtypes,params.xtime>0.2e6 & params.xtime<1e6),4),3);
% sortval                 = sortval(:,1)  - sortval(:,2);
% % [~,sortidx]             = sort(sortval);

% temp6                   = rand(size(temp))*1e-6;
% temp(temp==0)           = temp6(temp==0);
% temp(temp<=0)           = temp6(temp<=0);
% temp5                   = abs(temp5);

temp5(:,1)              = nanmean(nanmean(temp(idx,4,idx_ttypes_vis,params.xtime>0 & params.xtime<0.3e6),4),3);
temp5(:,2)              = nanmean(nanmean(temp(idx,5,idx_ttypes_aud,params.xtime>0e6 & params.xtime<0.2e6),4),3);
temp5(:,3)              = nanmean(nanmean(temp(idx,6,[idx_ttypes_vis idx_ttypes_aud],params.xtime>0.2e6 & params.xtime<1e6),4),3);

temp5                   = abs(temp5);
temp5(:,2)              = temp5(:,2)*50;
temp5(:,2)              = temp5(:,2).^2;
temp5(temp5(:,2)<1,2) = 1;

sortval                 = (-temp5(:,3) + temp5(:,1)) ./ temp5(:,2);
% sortval                 = temp5(:,1)  - temp5(:,3);

[~,sortidx]             = sort(sortval);

heatmat_vis             = heatmat_vis(sortidx,:,:);

%Make heatmap for auditory: trials::
heatmat_aud             = ones(sum(idx),params.nTimebins,3); %init
temp3                   = nanmean(temp(idx,5,idx_ttypes_aud,:),3);
heatmat_aud(:,:,2)      = heatmat_aud(:,:,2) - squeeze(temp3) / prctile(temp3(:),96);
heatmat_aud(:,:,3)      = heatmat_aud(:,:,3) - squeeze(temp3) / prctile(temp3(:),96);
temp4                   = nanmean(temp(idx,6,idx_ttypes_aud,:),3);
heatmat_aud(:,:,1)      = heatmat_aud(:,:,1) - squeeze(temp4) / prctile(temp4(:),96);
heatmat_aud(:,:,3)      = heatmat_aud(:,:,3) - squeeze(temp4) / prctile(temp4(:),96);

% heatmat_aud             = zeros(sum(idx),params.nTimebins,3);
% heatmat_aud(:,:,1)      = squeeze(nanmean(temp(idx,5,idx_trialtypes,:),3)) / max(nanmean(temp(idx,5,idx_trialtypes,:),3),[],'all');
% heatmat_aud(:,:,2)      = squeeze(nanmean(temp(idx,6,idx_trialtypes,:),3)) / max(nanmean(temp(idx,6,idx_trialtypes,:),3),[],'all');

%sort:
sortval(:,1)            = nanmean(nanmean(temp(idx,5,idx_ttypes_aud,params.xtime>0 & params.xtime<0.2e6),4),3);
sortval(:,2)            = nanmean(nanmean(temp(idx,6,idx_ttypes_aud,params.xtime>0.2e6 & params.xtime<1e6),4),3);
sortval                 = temp5(:,2)  - temp5(:,3);
% [~,sortidx]             = sort(sortval);
heatmat_aud             = heatmat_aud(sortidx,:,:);

figure; set(gcf,'units','normalized','Position',[0.05 0.59 0.33 0.3],'color','w'); hold all;
subplot(1,2,1); hold all;
imagesc(params.xtime,1:nNeurons,heatmat_vis)
xlim([params.t_pre params.t_post])
ylim([1 nNeurons])
plot([0 0],[1 nNeurons],'k:','LineWidth',0.5)
set(gca,'YTick',nNeurons,'XTick',-2e6:0.2e6:2e6,'XTickLabels',(-2e6:0.2e6:2e6)*1e-3)

subplot(1,2,2); hold all;
imagesc(params.xtime,1:nNeurons,heatmat_aud)
xlim([params.t_pre params.t_post])
ylim([1 nNeurons])
plot([0 0],[1 nNeurons],'k:','LineWidth',0.5)
set(gca,'YTick',[],'XTick',-2e6:0.2e6:2e6,'XTickLabels',(-2e6:0.2e6:2e6)*1e-3)

% export_fig(fullfile(params.savedir,sprintf('Zmat_overTime_Modelsplits')),'-eps','-nocrop')


%%






%% Depthcorrections:
tempcomb = {'20212019070104' 150; '20232019070170' -100; '20282019121322' -100; '20282019121626' -100; '20292019121240' -50; '20292019121343' -50};

for iSes = 1:size(tempcomb,1)
    spikeData.ChannelY(strcmp(spikeData.session_ID,tempcomb{iSes,1})) = spikeData.ChannelY(strcmp(spikeData.session_ID,tempcomb{iSes,1})) + tempcomb{iSes,2};
end

%%
params.depthcorrection          = -50;

spikeData.ChannelY              = spikeData.ChannelY + params.depthcorrection;

%% Make SUA heatmap over laminar depth for each trial type and for each model split:
%Histogram with binning on depth:

params.chdepthmin       = -12.5;
params.chdepthmax       = 1012.5;
params.resolution       = 25;

% binedges                = 0:25:1025;
binedges                = params.chdepthmin:params.resolution:params.chdepthmax;
nBins                   = length(binedges)-1;

params.finalYaxis       = binedges(1:end-1)+params.resolution/2;

nModelSplits            = size(ratemat,2);
SUAdepthmap_all         = zeros(nModelSplits,params.nSplits,nBins,params.nTimebins);

% iExp = 1;
%     idx_exp                 = ismember(spikeData.session_ID,sessionData.session_ID(strcmp(sessionData.Experiment,params.Experiments(iExp))));
idx_area                = strcmp(spikeData.area,'V1');% & idx_exp;

ratemat_z               = NaN(size(ratemat));
% ratemat_z               = NaN(size(ratemat));
for iNeuron = 1:nNeurons
    bsl = ratemat(iNeuron,:,:,params.xtime<0);
    ratemat_z(iNeuron,:,:,:) = (ratemat(iNeuron,:,:,:) - repmat(mean(bsl,4),1,1,1,params.nTimebins));
    ratemat_z(iNeuron,:,:,:) = ratemat_z(iNeuron,:,:,:) ./ repmat(nanstd(bsl(:)),size(ratemat_z(iNeuron,:,:,:)));
end

% ratemat_z               = ratemat;

% ratemat_z               = NaN(size(ratemat));
% for iNeuron = 1:nNeurons
%     bsl = ratemat(iNeuron,1,:,params.xtime<0);
%     ratemat_z(iNeuron,:,:,:) = (ratemat(iNeuron,:,:,:) - repmat(mean(bsl,4),1,6,1,params.nTimebins));
%     ratemat_z(iNeuron,:,:,:) = ratemat_z(iNeuron,:,:,:) ./ repmat(nanstd(bsl(:)),size(ratemat_z(iNeuron,:,:,:)));
% end


for iMod = 1:nModelSplits
    for iSplit = 1:params.nSplits
        for iBin = 1:nBins
            idx = idx_area & spikeData.ChannelY>=binedges(iBin) & spikeData.ChannelY<binedges(iBin+1);
            
            SUAdepthmap_all(iMod,iSplit,iBin,:) = nanmean(ratemat_z(idx,iMod,iSplit,:),1);
        end
    end
end

SUAdepthmap_all(isnan(SUAdepthmap_all)) = 0;

%Filter temporally and spatially:
% spat_filter     = fspecial('gaussian',[10 20],1.6);
% spat_filter     = fspecial('average',[10 2]);
spat_filter     = fspecial('average',[10 3]);

for iMod = 1:nModelSplits
    for iSplit = 1:params.nSplits
        temp = squeeze(SUAdepthmap_all(iMod,iSplit,:,:));
        SUAdepthmap_all(iMod,iSplit,:,:) = conv2(temp,spat_filter,'same');
    end
end

%%

plotModels      = [1 4 6 1 5 6];
plotTrialTypes  = {[1 2] [1 2] [1 2] [3 4] [3 4] [3 4]};
% plotTrialTypes  = {[1 2 3 4] [1 2 3 4] [1 2 3 4] [5 6 7 8] [5 6 7 8] [5 6 7 8]};

figure; set(gcf,'units','normalized','Position',[0.05 0.5 0.4 0.3],'color','w');
timeticks = [-0.2e6 0 0.2e6 0.4e6 0.6e6];
binticks = [0 250 500 750 1000];

cmapdata = getPyPlot_cMap('rainbow');
cmapdata = parula(128);

for iMod = 1:length(plotModels)
%     subplot(1,6,iMod);        hold all;
    subplot(2,3,iMod);        hold all;
    
    title(params.labels_cats(plotModels(iMod)),'FontSize',10);
    
    SUAdepthmap = squeeze(nanmean(SUAdepthmap_all(plotModels(iMod),plotTrialTypes{iMod},:,:),2));
    
    imagesc(params.xtime,params.finalYaxis,SUAdepthmap); 
    plot([0 0],[-2000 2000],':','Color',[1 1 1],'LineWidth',1)
    temp = caxis;
    caxis([-0.3 temp(2)])
    colormap(cmapdata)

    %Figure make up:
    if iMod==1 || iMod==4
        ylabel('Cortical depth (\mum)','FontSize', 9)
    end
    if iMod==5
        xlabel('Time (ms)','FontSize',9)
    end
    set(gca,'YDir','reverse') %correct for imagesc plotting inverse
    set(gca,'XTick',timeticks,'XTickLabel',timeticks*1e-3)
    set(gca,'YTick',binticks,'YTickLabel',binticks,'fontsize',7,'tickdir','out')
    xlim([find(params.xtime>-0.2e6,1) find(params.xtime<0.6e6,1,'last')])
    xlim([-0.2e6 0.6e6]);
    ylim([0 1000])
    
    if iMod==3
    temppos = get(gca,'position');
    colorbar('Location','eastoutside');
    set(gca,'Position',temppos)
%         cb=colorbar;
%         cb.Position = cb.Position + 1e-10;
    end
end

export_fig(fullfile(params.savedir,sprintf('LaminarHeatmap_ModelSplits_AV_3cohorts')),'-eps','-nocrop')



%% Video EV as a function of video PC:

%% Compute variance explained over all individual trials:
idx_time            = params.xtime>0 & params.xtime<1e6;
% idx_time            = params.xtime>0 & params.xtime<=0.2e6;

cvR2_full_videoPC       = NaN(nNeurons,params.nSVDs);
cvR2_full_videoPC_un    = NaN(nNeurons,params.nSVDs);

fprintf('Computing cvR2 for neuron        \n');

lastsesid = [];

startidx = find(strcmp(output.x_label,'videoMotion'),1);
startidx = startidx-1;

params.lambdastring         = 'lambda_1se';

for iNeuron = 1:nNeurons %loop over neurons
    fprintf(repmat('\b', 1, numel([num2str(iNeuron-1) num2str(nNeurons)])+2));
    fprintf('%d/%d\n',iNeuron,nNeurons);
    if ~isempty(output.modelFits(iNeuron,1).lambda)
        nTotalSpikeBins     = sum(~isnan(output.y(iNeuron,:))); %compute how many total spike bins there are from output struct
        nTrials             = nTotalSpikeBins / params.nTimebins; %how many trials iin this session
        %         idx_ses             = strcmp(sessionData.session_ID,spikeData.session_ID(iNeuron)); %get which session this is out of all ses predictors
        %         temptrialData       = MOL_getTempPerSes(sessionData.session_ID(idx_ses),trialData); %filter trialdata for this session
        temptrialData       = MOL_getTempPerSes(spikeData.session_ID(iNeuron),trialData);
        
        if ~(strcmp(spikeData.session_ID(iNeuron),lastsesid))
            load(fullfile(folderpath,sprintf('X_Full_GLM_%s.mat',spikeData.session_ID{iNeuron})))
            %  X_full              = squeeze(x(1,1:nTotalSpikeBins,1:params.nPredictors));
            lastsesid           = spikeData.session_ID(iNeuron);
        end
        
        %         X_full              = squeeze(output.x(idx_ses,1:nTotalSpikeBins,1:params.nPredictors)); %get predictor matrix for this session
        Y_full              = squeeze(output.y(iNeuron,1:nTotalSpikeBins))'; %get spike rates
        idx_time_all        = repmat(idx_time,1,nTrials); %make an index of all time bins that are included for computation of cvR2
        
%         Yh_full             = cvglmnetPredict(output.modelFits(iNeuron,1),X_full,params.lambdastring,'response'); %cv predicted firing rates for full model
%         temp_full           = 1 - var(Y_full(idx_time_all) - Yh_full(idx_time_all)) / var(Y_full(idx_time_all)); %compute R2
        
        for iSVD = 1:params.nSVDs
            idx                     = [(1:iSVD)+startidx (1:iSVD)+startidx+params.nSVDs ...
                                        (1:iSVD)+startidx+params.nSVDs*2];
            idx                     = ismember(1:params.nPredictors,idx);
            X_temp                  = X_full;
            X_temp(:,~idx)          = 0;
            Yh_svd                  = cvglmnetPredict(output.modelFits(iNeuron,1),X_temp,params.lambdastring,'response');
%             cvR2_cat(iNeuron,iC)    = var(Yh_cat(idx_time_all))/var(Y_full(idx_time_all));
            cvR2_full_videoPC(iNeuron,iSVD)    = 1 - var(Y_full(idx_time_all) - Yh_svd(idx_time_all)) / var(Y_full(idx_time_all));
            
            %compute firing rate with circularly shuffled predictor and compute uniquely contributed R2
            %for each SVD:
            iC = 4;
            %compute firing rate with circularly shuffled predictor and compute uniquely contributed R2:
            Yh_svd_shuf             = cvglmnetPredict(output.shuffleFits(iNeuron,iC),X_temp,params.lambdastring,'response');
            temp                    = 1 - var(Y_full(idx_time_all) - Yh_svd_shuf(idx_time_all)) / var(Y_full(idx_time_all));
            cvR2_full_videoPC_un(iNeuron,iSVD) =  cvR2_full_videoPC(iNeuron,iSVD) - temp;
        end
    end
end

%%

%%
% cvR2_full_videoPC2 = diff([ones(nNeurons,1) cvR2_full_videoPC],[],2);
cvR2_full_videoPC2 = diff(cvR2_full_videoPC,[],2);
cvR2_full_videoPC2 = cvR2_full_videoPC;
% cvR2_full_videoPC2_un = diff(cvR2_full_videoPC_un,[],2);
%%
% cvR2_full_videoPC2 = diff(cvR2_full_videoPC_un,[],2);

%% 
% cvR2_full_videoPC2(cvR2_full_videoPC2<prctile(cvR2_full_videoPC2(:),1)) = prctile(cvR2_full_videoPC2(:),1);
% cvR2_full_videoPC2(cvR2_full_videoPC2>prctile(cvR2_full_videoPC2(:),99.9)) = prctile(cvR2_full_videoPC2(:),99);

cvR2_full_videoPC2(cvR2_full_videoPC2<0) = 0;

%%
params.areas            = {'V1'};% 'PPC' 'CG1'};
params.nAreas           = length(params.areas); 

figure; set(gcf,'units','normalized','Position',[0.2 0.43 0.22 0.3],'color','w'); hold all;

params.colors_splits = {[0.4 0.4 0.4] [0 0 0.8] [0.8 0 0] [0.1 0.7 0.3]};
handles = [];
for iArea = 1:params.nAreas
    idx     = strcmp(spikeData.area,params.areas{iArea});
%     plot([0 mean(cvR2_full_videoPC2(idx,:),1)])
%      plot(cvR2_full_videoPC2(idx,:)','Color',[0.6 0.6 0.6],'LineWidth',0.25)
    h = errorbar(mean(cvR2_full_videoPC2(idx,:),1),std(cvR2_full_videoPC2(idx,:),[],1) / sqrt(nNeurons),'Color',[0 0.1 1],'LineWidth',2);
%     h = errorbar(mean(cvR2_full_videoPC2_un(idx,:),1),std(cvR2_full_videoPC2_un(idx,:),[],1) / sqrt(nNeurons),'Color',[0 0.1 1],'LineWidth',2);
    errorbar_tick(h,0.001,'units');
end
plot([0 500],[0.0003 0.0003],'k:','LineWidth',0.5)

ylabel('Explained Variance')
xlabel('Video PC #')
xlim([0.5 params.nSVDs])
ylim([0 0.01])
set(gca,'XTick',[1 5 10 15 20 25 30],'YTick',[0 0.0025 0.005 0.0075 0.01])

% export_fig(fullfile(params.savedir,sprintf('EV_vs_videoPC_V1neurons')),'-eps','-nocrop')

%% Control figure, volume vs cvR2 
params.colors_cats      = {[0.5 0.2 0.3] [0 0 0.8] [0.8 0 0] [0.1 0.7 0.3]};

spikeData.dblSoundLevel = NaN(nNeurons,1);
for iSes = 1:nSessions
    spikeData.dblSoundLevel(strcmp(spikeData.session_ID,sessionData.session_ID(iSes))) = sessionData.dblSoundLevel(iSes);
end

figure; set(gcf,'units','normalized','Position',[0.2 0.43 0.22 0.3],'color','w'); hold all;
subplot(2,1,1); hold all;
h = scatter(spikeData.dblSoundLevel+randn(size(spikeData.dblSoundLevel))*1e-3,cvR2_cat(:,3),15,params.colors_cats{3},'.');
xlim([0.5 0.85])
ylim([0 0.4])
ylabel('Aud. EV')
[r,p]               = corr(spikeData.dblSoundLevel,cvR2_cat(:,3),'rows','complete','type', 'Spearman');
fprintf('Spearman rank correlation between sound volume \n and auditory-related EV, r=%1.3f, p=%1.3f; \n',r,p)
subplot(2,1,2); hold all;
h = scatter(spikeData.dblSoundLevel+randn(size(spikeData.dblSoundLevel))*1e-3,cvR2_cat(:,4),15,params.colors_cats{4},'.');
xlim([0.5 0.85])
ylim([0 0.4])
ylabel('Motor EV')
xlabel('Session sound volume (dBA)')

[r,p]               = corr(spikeData.dblSoundLevel,cvR2_cat(:,4),'rows','complete','type', 'Spearman');
fprintf('Spearman rank correlation between sound volume \n and motor-related EV, r=%1.3f, p=%1.3f; \n',r,p)
export_fig(fullfile(params.savedir,sprintf('EV_vs_soundvolume_scatter_Aud_Motor')),'-eps','-nocrop')

%% 
% % %% Correlation between cvR2 of auditory, visual, and movement variables:
% % 
% % figure; set(gcf,'units','normalized','Position',[0.05 0.59 0.25 0.17],'color','w'); hold all;
% % subplot(1,2,1); hold all;
% % 
% % idx_exp                 = ismember(spikeData.session_ID,sessionData.session_ID(ismember(sessionData.Experiment,params.Experiments)));
% % 
% % idx_area                = strcmp(spikeData.area,'V1');
% % idx                     = idx_exp & idx_area;
% % 
% % xdata                   = cvR2_cat(idx,2);
% % ydata                   = cvR2_cat(idx,4);
% % 
% % scatter(xdata,ydata,30,'bo','filled','MarkerFaceAlpha',0.4) %Color',params.colors_cats{1},'LineWidth',0.5);
% % 
% % ylim([0 0.4]); xlim([0 0.4])
% % 
% % %statistics:
% % tbl             = table(ydata,xdata,G_mou(idx),'VariableNames',{'Y','X','Mouse'}); %Create table for mixed model
% % lme             = fitlme(tbl,'Y~X+(1|Mouse)'); %construct linear mixed effects model with fixed effect of temporal window (NO random intercept for mice as cohorts are tested)
% % stats           = dataset2table(anova(lme,'DFMethod','Satterthwaite')); %Perform ANOVA on model and output as matrix
% % r               = corr(xdata,ydata,'rows','complete');
% % 
% % fprintf('Variance explained was uncorrelated between visual and motor predictors \n')
% % fprintf('R=%1.2f, F(%d,%2.0f)=%1.2f, p=%1.2f; \n',r,stats{2,3},stats{2,4},stats{2,2},stats{2,5})
% % 
% % subplot(1,2,2); hold all;
% % 
% % xdata                   = cvR2_cat(idx,3);
% % ydata                   = cvR2_cat(idx,4);
% % 
% % %statistics:
% % tbl             = table(ydata,xdata,G_mou(idx),'VariableNames',{'Y','X','Mouse'}); %Create table for mixed model
% % lme             = fitlme(tbl,'Y~X+(1|Mouse)'); %construct linear mixed effects model with fixed effect of temporal window (NO random intercept for mice as cohorts are tested)
% % stats           = dataset2table(anova(lme,'DFMethod','Satterthwaite')); %Perform ANOVA on model and output as matrix
% % r               = corr(xdata,ydata,'rows','complete');
% % 
% % fprintf('and weakly correlated between auditory and motor predictors \n')
% % fprintf('r=%1.2f, F(%d,%2.0f)=%1.2f, p=%1.2e; \n',r,stats{2,3},stats{2,4},stats{2,2},stats{2,5})
% % 
% % scatter(xdata,ydata,30,'ro','filled','MarkerFaceAlpha',0.4) %Color',params.colors_cats{1},'LineWidth',0.5);
% % ylim([0 0.4]); xlim([0 0.1])
% % 
% % export_fig(fullfile(params.savedir,sprintf('Scatter_cvR2_AV_3cohorts')),'-eps','-nocrop')



% %% Set as significant if responding to one the two stimuli:
% 
% signmat_any(:,1,:) = any(signmat(:,[1 2],:),2);
% signmat_any(:,2,:) = any(signmat(:,[3 4],:),2);
% signmat_any(:,3,:) = any(signmat(:,[5 6],:),2);
% signmat_any(:,4,:) = any(signmat(:,[7 8],:),2);
% 
% signmat_any = signmat;
% 
% %% Show fracion of significantly responsive neurons when regressing out movement variability:
% 
% params.colors_trialtypes    = [params.colors_visual_opto(1:2) params.colors_audio_opto(1:2)];
% params.colors_trialtypes    = [params.colors_visual_opto([2 1]) params.colors_audio_opto([2 1])];
% 
% params.colors_trialtypes    = [params.colors_ztrials(1:2) {[0.9 0.3 0.1] [0.5 0 0.1]}];
% 
% figure; set(gcf,'units','normalized','Position',[0.05 0.4 0.33 0.33],'color','w'); hold all;
% 
% % iExp = 3;
% %     idx_exp                 = ismember(spikeData.session_ID,sessionData.session_ID(strcmp(sessionData.Experiment,params.Experiments(iExp))));
% idx_exp                 = ismember(spikeData.session_ID,sessionData.session_ID(ismember(sessionData.Experiment,params.Experiments)));
% 
% 
% idx_area    = strcmp(spikeData.area,'V1');
% % idx_area    = strcmp(spikeData.area,'V1') & idx_exp;
% 
% idx_models  = [1 4];
% 
% labels      =  {'Raw' '-Motor'};
% 
% colors_models = {[0 0 0] [0.4 0.6 0.8] [0.3 0.4 0.1] [0.8 0.2 0.6]};
% 
% subplot(1,2,1); hold all;
% for iSplit = 1:4
%     datatoplot = squeeze(sum(signmat_any(idx_area,iSplit,idx_models),1)) / sum(idx_area);
%     h = bar((1:length(idx_models))+(iSplit-1)*2,datatoplot,'FaceColor',params.colors_trialtypes{iSplit});
% end
% set(gca,'XTick',1:8,'XTickLabels',repmat(labels,1,4),'YTick',0:0.2:1);
% ylim([0 0.6])
% title('V1')
% 
% idx_area    = strcmp(spikeData.area,'A1');
% subplot(1,2,2); hold all;
% for iSplit = 1:4
%     datatoplot = squeeze(sum(signmat_any(idx_area,iSplit,idx_models),1)) / sum(idx_area);
%     h = bar((1:length(idx_models))+(iSplit-1)*2,datatoplot,'FaceColor',params.colors_trialtypes{iSplit});
% end
% set(gca,'XTick',1:8,'XTickLabels',repmat(labels,1,4),'YTick',0:0.2:1);
% ylim([0 0.8])
% legend(params.labels_splits)
% title('AC')
% 
% %% 
% 



