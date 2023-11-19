
%%
startover

%% Load dataset:
savedate                    = '09-12-21';
% savedate                    = '03-01-23cross';
savedate                    = '20-12-22';
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
% output.shuffleFits  = loadstruct.output.shuffleFits;
output.modelFits    = loadstruct.output.modelFits;
output.x_label      = loadstruct.output.x_label;
nTotalneurons       = 0;

sessionData         = struct();
trialData           = struct();
spikeData           = struct();
% videoData           = struct();

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
%     videoData           = AppendStruct(videoData,loadstruct.videoData);
    
    nNeurons = length(loadstruct.spikeData.session_ID);
    idx = nTotalneurons+1:nTotalneurons+nNeurons;
    
%     output.x(iF,:,:)    = loadstruct.output.x;
    output.x_sesid(iF,1)= loadstruct.output.x_sesid;
    output.y(idx,:)     = loadstruct.output.y;
    if iF>1
        output.modelFits    = [output.modelFits; loadstruct.output.modelFits];
        output.shuffleFits  = [output.shuffleFits; loadstruct.output.shuffleFits];
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

%% %remove animal 1012 with recordings in LM according to histology
spikeData.area(strcmp(spikeData.session_ID,'10122019041022') & strcmp(spikeData.area,'A1')) = deal({'LM'});
spikeData.area(strcmp(spikeData.session_ID,'10122019041126') & strcmp(spikeData.area,'A1')) = deal({'LM'});

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

%% FIGURES: 

%% Rename splits and compute for pupil separately:
params.shuffleVars{4} = {'videoMotion'};
params.shuffleVars{5} = {'pupilArea' 'pupilX'   'pupilY'};

params.nShuffleCats = 5;

%% Compute variance explained over all individual trials:
idx_time            = params.xtime>0 & params.xtime<=0.2e6;

cvR2_full           = NaN(nNeurons,1);
cvR2_cat            = NaN(nNeurons,params.nShuffleCats);
% cvR2_cat_un         = NaN(nNeurons,params.nShuffleCats);

fprintf('Computing cvR2 for neuron        \n');

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
            lastsesid           = spikeData.session_ID(iNeuron);
        end
        
        Y_full              = squeeze(output.y(iNeuron,1:nTotalSpikeBins))'; %get spike rates
        Yh_full             = cvglmnetPredict(output.modelFits(iNeuron,1),X_full,params.lambdastring,'response'); %cv predicted firing rates for full model
        idx_time_all        = repmat(idx_time,1,nTrials); %make an index of all time bins that are included for computation of cvR2
        cvR2_full(iNeuron,1) = 1 - var(Y_full(idx_time_all) - Yh_full(idx_time_all)) / var(Y_full(idx_time_all)); %compute R2

        for iC = 1:params.nShuffleCats
            %Predict firing rate when using only one predictor set, compute R2 and store value:
            idx                     = ismember(output.x_label,params.shuffleVars{iC});
            X_temp                  = X_full;
            X_temp(:,~idx)          = 0;
            Yh_cat                  = cvglmnetPredict(output.modelFits(iNeuron,1),X_temp,params.lambdastring,'response');
            cvR2_cat(iNeuron,iC)    = 1 - var(Y_full(idx_time_all) - Yh_cat(idx_time_all)) / var(Y_full(idx_time_all));
        end
    end
end

%% Show cvR2 for each predictor subset and area:
params.nExperiments     = length(params.Experiments);
params.areas            = {'V1' 'A1'};% 'PPC' 'CG1'};
params.nAreas           = length(params.areas); 

figure; set(gcf,'units','normalized','Position',[0.2 0.43 0.22 0.3],'color','w'); hold all;

params.colors_splits = {[0.4 0.4 0.4] [0 0 0.8] [0.8 0 0] [0.1 0.7 0.3]};
params.colors_splits = {[0.4 0.4 0.4] [0 0 0.8] [0.8 0 0] [0.1 0.7 0.3] [0.7 0.2 0.8]};
params.Yvarlabels    = {'True' 'Full model' 'Trial' 'Visual' 'Audio' 'Motor' 'Pupil'};

handles = [];

for iArea = 1:params.nAreas
    idx     = strcmp(spikeData.area,params.areas{iArea});
    
    for iC = 1:params.nShuffleCats
        tmp     = cvR2_cat(idx,iC);
%         tmp     = cvR2_cat_un(idx,iC);
        handles(iC) = bar(iArea + iC/(params.nShuffleCats+1),nanmean(tmp),0.18,'k');
        set(handles(iC),'FaceColor',params.colors_splits{iC})
        errorbar(iArea + iC/(params.nShuffleCats+1),nanmean(tmp),nanstd(tmp)/sqrt(sum(idx)),'k','LineWidth',1,'CapSize',6);
    end
end

tmp     = cvR2_cat(strcmp(spikeData.area,'V1'),5);
fprintf('Pupil size and location explained only a minor fraction of V1 variance (%1.4f +- %1.2e, mean +- sem) and was therefore not analyzed separately\n',nanmean(tmp),nanstd(tmp)/sqrt(sum(idx)))

%statistics:
% Y               = cvR2_cat_un(:);
Y               = cvR2_cat(:);
X_var           = repmat(params.Yvarlabels(3:7),nNeurons,1);
X_var           = X_var(:);
X_area          = repmat(spikeData.area,1,5);
X_area          = X_area(:);
G_mou2          = repmat(G_mou,1,5);
G_mou2          = G_mou2(:);

idx_V1          = strcmp(X_area,'V1');
idx_AC          = strcmp(X_area,'A1');

tbl             = table(Y(idx_V1),X_var(idx_V1),G_mou2(idx_V1),'VariableNames',{'Y','Var','Mouse'}); %Create table for mixed model
lme             = fitlme(tbl,'Y~Var+(1|Mouse)'); %construct linear mixed effects model with fixed effect of temporal window (NO random intercept for mice as cohorts are tested)
stats           = dataset2table(anova(lme,'DFMethod','Satterthwaite')); %Perform ANOVA on model and output as matrix
fprintf('F(%d,%2.0f)=%1.2f, p=%1.2e; \n',stats{2,3},stats{2,4},stats{2,2},stats{2,5})

for iC = 1:5
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
legend(handles(1:5),{'Trial' 'Visual' 'Audio' 'Motor' 'Pupil'},'Location','NorthWest'); legend boxoff;

% export_fig(fullfile(params.savedir,sprintf('Bar_cvR2_V1_AC')),'-eps','-nocrop')
