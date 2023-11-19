%% Oude Lohuis et al. 2023 Nat Neurosci
% Triple dissociation of auditory, visual, and motor processing in primary visual cortex
% MOL (C) 2023

startover

%% Parameter settings:
params                      = params_histresponse_auV1;% All time is in microseconds

params.Experiments          = {'ChangeDetectionConflictDecor' 'VisOnlyTwolevels' 'ChangeDetectionConflict' }; %Which versions of the task to load data from
params.ExperimentLabels     = {'NE' 'UST' 'MST'}; %Which versions of the task to load data from
params.AlignOn              = 'stimChange';         %On which timestamp to align trials to
params.nExperiments         = length(params.Experiments);
params                      = MOL_getColors_CHDET(params);

params.videofield           = 'motSVD';

params.smoothSVD            = 0;
params.nSVDs                = 25;

params.areas                = {'V1'};
params.nAreas               = length(params.areas);

%Parameters for AUC:
params.minTrialCond         = 10;              %Number of minimum trials that both conditions need to have to be discriminated
params.smoothing            = 0;
params.zscore               = 0;
params.subtr_baseline       = 1;               %Subtract baseline or not
params.nshuffle             = 1000;            %Number of shuffles to base permutation test on
params.alpha                = 0.05;            %Significance level for permutation test

params.savedir              = 'E:\Documents\PhD\Figures\Project CHDET\Results - auV1\6FeatureCoding'; %output saving dir for figures

%% Get data:
[Data] = MOL_GetData('E:','CHDET',params.Experiments,{},[],{'sessionData' 'trialData_newtrials' 'spikeData'  'videoData'});
sessionData         = Data.sessionData;
trialData           = Data.trialData;
spikeData           = Data.spikeData;
videoData           = Data.videoData;

%% Remove last 20 trials:
trialData           = MOL_RemoveLastnTrials(trialData,20);

%% Filter out neurons based on quality:
spikeData           = MOL_filterNeurons(sessionData,trialData,spikeData);

%% Remove sessions that are passive: 
sesids                      = sessionData.session_ID(strcmp(sessionData.State,'Behaving'));
fprintf('Removed %d/%d sessions with passive stimulation\n',length(sessionData.session_ID)-length(sesids),length(sessionData.session_ID));
[sessionData, trialData, spikeData,videoData]        = MOL_getTempPerSes(sesids,sessionData,trialData,spikeData,videoData);

%% Filter out sessions with no systematic post orientation or frequency:
sesids = unique(trialData.session_ID(~isnan(trialData.visualOriPostChangeNorm) & ~isnan(trialData.audioFreqPostChangeNorm)));
[sessionData, trialData, spikeData,videoData]        = MOL_getTempPerSes(sesids,sessionData,trialData,spikeData,videoData);

%% Filter out sessions with muscimol:
sesids              = sessionData.session_ID(cellfun(@isempty,sessionData.MuscimolArea));
fprintf('Removed %d/%d sessions  with muscimol manipulations\n',length(sessionData.session_ID)-length(sesids),length(sessionData.session_ID));
[sessionData, trialData, spikeData,videoData]        = MOL_getTempPerSes(sesids,sessionData,trialData,spikeData,videoData);

%% Filter out trials with photostimulation in V1:
sesids              = sessionData.session_ID(sessionData.UseOpto & strcmp(sessionData.PhotostimArea,'V1'));
idx                 = trialData.hasphotostim==1 & ismember(trialData.session_ID,sesids);
trialFields         = fieldnames(trialData);
fprintf('Removed %d/%d trials  with photostim in V1\n',sum(idx),length(trialData.session_ID));
for iF = 1:length(trialFields)
    trialData.(trialFields{iF}) = trialData.(trialFields{iF})(~idx);
end

%% Filter neurons on area:
spikeFields     = fieldnames(spikeData);
idx             = ismember(spikeData.area,params.areas);
fprintf('Filtered %d/%d neurons based on area\n',sum(idx),length(spikeData.session_ID));
for iF = 1:length(spikeFields)
    spikeData.(spikeFields{iF}) = spikeData.(spikeFields{iF})(idx);
end

%% Get only sessions with neurons recorded:
[sessionData,trialData,spikeData,videoData] = MOL_getTempPerSes(unique(spikeData.session_ID),sessionData,trialData,spikeData,videoData);

%% Filter sessions based on containing motSVD data:
% [videoData] = MOL_getTempPerSes(videoData.session_ID(~cellfun(@isempty,videoData.motSVD)),videoData);
[sessionData,trialData,spikeData,videoData] = MOL_getTempPerSes(videoData.session_ID(~cellfun(@isempty,videoData.motSVD)),sessionData,trialData,spikeData,videoData);

%% Save dataset:
save('Dataset2_1.mat','params','sessionData','trialData','videoData','-v7.3')

%% Or load dataset:
load Dataset2_1.mat

%% Report dataset:
fprintf('Dataset: %d sessions, %d trials, %d neurons, %d videos\n',length(sessionData.session_ID),length(trialData.session_ID),length(spikeData.session_ID),length(videoData.session_ID));
nSessions   = length(sessionData.session_ID);
nNeurons    = length(spikeData.session_ID);

%% 






%% Make psth of mot svd from video:
[hist_mat,hist_mat_tot,hist_mat_z,params] = MOL_PSTH_video(sessionData,trialData,videoData,params);

nTrials             = length(trialData.stimChange);

%% Compute AUC for neurons:
[AUC_ORI_CELL,AUC_FREQ_CELL,pVal_ORI_CELL,pVal_FREQ_CELL,respVis,respAud] = calc_AUC_spikes(params,sessionData,trialData,spikeData);

%% Compute AUC for video ME:
[AUC_ORI_VID,AUC_FREQ_VID,pVal_ORI_VID,pVal_FREQ_VID] = calc_AUC_video(params,trialData,videoData,hist_mat_z);

%% Make figure of two example sessions where responses to the tones are correlated within session for frequency, but not orientation:
close all;
params.zscore = 1;

params.exampleSes = {'20452021042245' '20112018080800'};

% for iSes = 1:nSessions
for iSes = find(ismember(sessionData.session_ID,params.exampleSes))'
    idx             = strcmp(spikeData.session_ID,sessionData.session_ID(iSes));
    nSesNeurons     = sum(idx);
    if nSesNeurons>15
        tempAUCORI      = AUC_ORI_CELL(idx);
        tempAUCFRQ      = AUC_FREQ_CELL(idx);
        
        %Get the relevant data for each session individually:
        [temptrialData,tempspikeData]        = MOL_getTempPerSes(sessionData.session_ID(iSes),trialData,spikeData);
        %Get the right indices:
        splits                  = {};
        splits{1}       = strcmp(temptrialData.trialType,'X') & temptrialData.visualOriChangeNorm==3 & ismember(temptrialData.visualOriPostChangeNorm,[1 2]);
        splits{2}       = strcmp(temptrialData.trialType,'X') & temptrialData.visualOriChangeNorm==3 & ismember(temptrialData.visualOriPostChangeNorm,[3 4]);
        splits{3}       = strcmp(temptrialData.trialType,'Y') & temptrialData.audioFreqChangeNorm==3 & ismember(temptrialData.audioFreqPostChangeNorm,[1 2]);
        splits{4}       = strcmp(temptrialData.trialType,'Y') & temptrialData.audioFreqChangeNorm==3 & ismember(temptrialData.audioFreqPostChangeNorm,[3 4]);
        
        tempresp        = NaN(nSesNeurons,4);
        %Compute histogram:
        for iNeuron = 1:nSesNeurons
            events_ts               = temptrialData.(params.AlignOn);
            hist_mat                = calc_psth(events_ts,tempspikeData.ts{iNeuron},params);    %Construct histogram matrix:
            for iSplit = 1:4
                tempresp(iNeuron,iSplit)    = nanmean(nanmean(hist_mat(splits{iSplit},params.xtime>params.twin_resp_start & params.xtime<params.twin_resp_stop),2),1);
            end
        end
        
        [~,sortidx] = sort(abs(diff(tempresp(:,[1 2]),[],2)));
        tempresp(:,[1 2]) = tempresp(sortidx,[1 2]);
        [~,sortidx] = sort(abs(diff(tempresp(:,[3 4]),[],2)));
        tempresp(:,[3 4]) = tempresp(sortidx,[3 4]);
        
        figure; set(gcf,'units','normalized','Position',[0.25 0.47 0.14 0.3],'color','w'); hold all
        subplot(1,2,1); hold all
        imagesc(tempresp(:,[1 2]));
        caxis([-0.3 1.4])
        ylim([0.5 nSesNeurons+0.5])
        set(gca,'XTick',[],'YTick',[])
        title(sessionData.session_ID{iSes})
        colorbar()

        subplot(1,2,2); hold all
        imagesc(tempresp(:,[3 4]))
        caxis([-0.3 1.4])
        ylim([0.5 nSesNeurons+0.5])
        set(gca,'XTick',[],'YTick',[])
        colorbar()
        
    end
end

%% For multi-level statistics:
G_mou           = cell(nNeurons,1);
uMice           = unique(sessionData.mousename);
for iMouse = 1:length(uMice)
    G_mou(ismember(spikeData.session_ID,sessionData.session_ID(strcmp(sessionData.mousename,uMice{iMouse})))) = uMice(iMouse);
end

%% Make figure of the correlation between neuronal and video AUC:

temp_AUC_ORI_vid_cells = NaN(size(spikeData.session_ID));
temp_AUC_FREQ_vid_cells = NaN(size(spikeData.session_ID));
for iSes = 1:nSessions
    temp_AUC_ORI_vid_cells(strcmp(spikeData.session_ID,sessionData.session_ID(iSes)))       = AUC_ORI_VID(iSes,1);
    temp_AUC_FREQ_vid_cells(strcmp(spikeData.session_ID,sessionData.session_ID(iSes)))      = AUC_FREQ_VID(iSes,1);
end

figure; set(gcf,'units','normalized','Position',[0.25 0.47 0.4 0.3],'color','w'); hold all
subplot(1,2,1); hold all;
scatter(AUC_ORI_CELL,temp_AUC_ORI_vid_cells,25,[0 0 0.8],'filled');
xlim([-1 1]); ylim([-1 1]);
plot([0 0],[-1 1],'k:','LineWidth',0.5)
plot([-1 1],[0 0],'k:','LineWidth',0.5)
set(gca,'XTick',[-1 0 1],'YTick',[-1 0 1]);
xlabel('AUC (neurons)'); ylabel('AUC (video ME)')
title('Orientation')

%statistics:
tbl             = table(AUC_ORI_CELL,temp_AUC_ORI_vid_cells,G_mou,'VariableNames',{'AUC_cell','AUC_vid','Mouse'}); %Create table for mixed model
lme             = fitlme(tbl,'AUC_cell~AUC_vid+(1|Mouse)'); %construct linear mixed effects model with fixed effect of temporal window (NO random intercept for mice as cohorts are tested)
stats           = dataset2table(anova(lme,'DFMethod','Satterthwaite')); %Perform ANOVA on model and output as matrix
fprintf('Correlation orientation selectivity spikes and video: \n')
fprintf('F(%d,%2.0f)=%1.2f, R=%1.3f, p=%1.2f; \n',stats{2,3},stats{2,4},stats{2,2},sqrt(lme.Rsquared.Ordinary),stats{2,5})
text(-0.9,0.7,sprintf('R=%1.3f\np=%1.3f',sqrt(lme.Rsquared.Ordinary),stats{2,5}))

subplot(1,2,2); hold all;
scatter(AUC_FREQ_CELL,temp_AUC_FREQ_vid_cells,25,[0.8 0 0],'filled');
xlim([-1 1]); ylim([-1 1]);
plot([0 0],[-1 1],'k:','LineWidth',0.5)
plot([-1 1],[0 0],'k:','LineWidth',0.5)
set(gca,'XTick',[-1 0 1],'YTick',[-1 0 1]);
xlabel('AUC (neurons)'); ylabel('AUC (video ME)')
title('Frequency')

%statistics:
tbl             = table(AUC_FREQ_CELL,temp_AUC_FREQ_vid_cells,G_mou,'VariableNames',{'AUC_cell','AUC_vid','Mouse'}); %Create table for mixed model
lme             = fitlme(tbl,'AUC_cell~AUC_vid+(1|Mouse)'); %construct linear mixed effects model with fixed effect of temporal window (NO random intercept for mice as cohorts are tested)
stats           = dataset2table(anova(lme,'DFMethod','Satterthwaite')); %Perform ANOVA on model and output as matrix
fprintf('Correlation frequency selectivity spikes and video: \n')
fprintf('F(%d,%2.0f)=%1.2f, R=%1.3f, p=%1.2e; \n',stats{2,3},stats{2,4},stats{2,2},sqrt(lme.Rsquared.Ordinary),stats{2,5})
text(-0.9,0.7,sprintf('R=%1.3f\np=%1.3e',sqrt(lme.Rsquared.Ordinary),stats{2,5}))

export_fig(fullfile(params.savedir,sprintf('Scatter_AUC_cell_vs_video')),'-eps','-nocrop')

