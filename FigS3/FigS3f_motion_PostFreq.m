%% Oude Lohuis et al. 2023 Nat Neurosci
% Triple dissociation of auditory, visual, and motor processing in primary visual cortex
% MOL (C) 2023

%Script shows poststimulus video motion for two different stimuli

startover

%% Parameter settings
params                      = params_histresponse_auV1;% All time is in microseconds

params.Experiments          = {'ChangeDetectionConflictDecor' 'VisOnlyTwolevels' 'ChangeDetectionConflict' }; %Which versions of the task to load data from
params.ExperimentLabels     = {'NE' 'UST' 'MST'}; %Labels for the different experiments
params.nExperiments         = length(params.Experiments);

params.AlignOn              = 'stimChange';      %On which timestamp to align as t=0

params.videofield           = 'motSVD';

params.area                 = 'V1';

params.fs                   = 25; %Hz

params.t_pre                = -1e6;
params.t_post               = 4e6;

params                      = MOL_getColors_CHDET(params);

params.savedir              = 'E:\Documents\PhD\Figures\Project CHDET\Results - auV1\6FeatureCoding\MovementTuning';

%% Get input arguments:
[Data] = MOL_GetData('E:','CHDET',params.Experiments,{'2009' '2012' '2013' '2011'},{},{'sessionData' 'trialData_newtrials' 'videoData'});
sessionData     = Data.sessionData;
trialData       = Data.trialData;
videoData       = Data.videoData;

% Remove last 20 trials:
trialData = MOL_RemoveLastnTrials(trialData,20);

%% Filter out trials with photostimulation in V1:
sesids              = sessionData.session_ID(sessionData.UseOpto & strcmp(sessionData.PhotostimArea,'V1'));
idx                 = trialData.hasphotostim==1 & ismember(trialData.session_ID,sesids);
trialFields         = fieldnames(trialData);
fprintf('Removed %d/%d trials  with photostim in V1\n',sum(idx),length(trialData.session_ID));
for iF = 1:length(trialFields)
    trialData.(trialFields{iF}) = trialData.(trialFields{iF})(~idx);
end

%% Filter out sessions with no systematic post orientation or frequency:
sesids = unique(trialData.session_ID(~isnan(trialData.visualOriPostChangeNorm) & ~isnan(trialData.audioFreqPostChangeNorm)));
[sessionData, trialData, videoData]        = MOL_getTempPerSes(sesids,sessionData,trialData,videoData);

%% Only sessions with video:
% Filter sessions based on containing videoData:
[sessionData,trialData,videoData] = MOL_getTempPerSes(unique(videoData.session_ID),sessionData,trialData,videoData);

% Filter sessions based on containing motSVD data:
[sessionData,trialData,videoData] = MOL_getTempPerSes(videoData.session_ID(~cellfun(@isempty,videoData.motSVD)),sessionData,trialData,videoData);

%% Report dataset:
fprintf('Dataset: %d sessions, %d trials, %d videos \n',length(sessionData.session_ID),length(trialData.session_ID),length(videoData.session_ID));
nSessions = length(sessionData.session_ID);

%%



%% Set variables for video processing:
params.videofield       = 'motSVD';

edges                   = params.t_pre:1e6/params.fs:params.t_post;
params.xtime_video      = edges+1e6/params.fs/2;
params.nTimebins_video  = numel(edges);
nTrials             	= length(trialData.stimChange);

%% Make figure of mean motion response for the different orientations and frequencies:

params.exSes            = {'20092018082260'}; %Example with similar movements to different stimuli 
params.exSes            = {'20112018081007'}; %Example with different movements to different frequencies
% params.exSes            = {'20132018082340'}; %example with different movements to different orientations

params.colors_splits    = {[255 147 12] [31 119 180] [255 147 12] [31 119 180]};
params.colors_splits       = cellfun(@(x) x/256,params.colors_splits,'UniformOutput',false);

params.nSVDs            = 3; %amount of pca's to plot

params.nSplits          = 4;

params.transparent      = 0;

params.triallabels      = {'V_A_B' 'V_C_D' 'A_A_B' 'A_C_D'};

meantoplot      = NaN(params.nTimebins_video,params.nSVDs,params.nSplits); %Init matrix vars
errortoplot     = NaN(params.nTimebins_video,params.nSVDs,params.nSplits);

[tempsessionData,temptrialData,tempvideoData] = MOL_getTempPerSes(params.exSes,sessionData,trialData,videoData);
splits          = {};
splits{1}       = strcmp(temptrialData.trialType,'X') & temptrialData.visualOriChangeNorm==3 & ismember(temptrialData.visualOriPostChangeNorm,[1 2]);
splits{2}       = strcmp(temptrialData.trialType,'X') & temptrialData.visualOriChangeNorm==3 & ismember(temptrialData.visualOriPostChangeNorm,[3 4]);
splits{3}       = strcmp(temptrialData.trialType,'Y') & temptrialData.audioFreqChangeNorm==3 & ismember(temptrialData.audioFreqPostChangeNorm,[1 2]);
splits{4}       = strcmp(temptrialData.trialType,'Y') & temptrialData.audioFreqChangeNorm==3 & ismember(temptrialData.audioFreqPostChangeNorm,[3 4]);

nTempTrials         = length(temptrialData.stimChange);

hist_mat_exses      = NaN(params.nTimebins_video,params.nSVDs,nTempTrials);

fprintf('Computing svd response for trial        \n');
for iTrial = 1:nTempTrials
    fprintf(repmat('\b', 1, numel([num2str(iTrial-1) num2str(nTempTrials)])+1));
    fprintf('%d/%d',iTrial,nTempTrials);
    hist_mat_exses(:,:,iTrial)                  = interp1(tempvideoData.ts{1}-temptrialData.(params.AlignOn)(iTrial),tempvideoData.(params.videofield){1}(:,1:params.nSVDs),params.xtime_video,'linear');
end
fprintf('\n')

for iSplit = 1:4
    meantoplot(:,:,iSplit)      = nanmean(hist_mat_exses(:,:,splits{iSplit}),3);
    errortoplot(:,:,iSplit)     = nanstd(hist_mat_exses(:,:,splits{iSplit}),[],3) / sqrt(sum(splits{iSplit}));
end

%Plot figure:
figure; hold all; set(gcf,'units','normalized','Position',[0.05 0.35 0.4 0.14],'color','w');

%visual first:
for iSVD = 1:params.nSVDs
    subplot(1,6,iSVD); hold all;
    handles = [];
    for iSplit = 1:2
%         handles(iSplit) = plot(params.xtime_video*1e-6,squeeze(meantoplot(:,iSVD,iSplit)),'Color',params.colors_splits{iSplit},'LineWidth',1);
        
        h = shadedErrorBar(params.xtime_video*1e-6,squeeze(meantoplot(:,iSVD,iSplit)),squeeze(errortoplot(:,iSVD,iSplit)),{'-','markerfacecolor',params.colors_splits{iSplit},'LineWidth',1},params.transparent);
        delete(h.edge(1:2));
        h.mainLine.Color = params.colors_splits{iSplit};    h.patch.FaceColor = params.colors_splits{iSplit};
        handles(iSplit) = h.mainLine; 
    end
    %Figure make up:
    if iSVD==2
        title('Visual')
    end
    curylim = get(gca,'ylim');
    plot([0 0],curylim,'-','Color',[0.6 0.6 0.6],'LineWidth',0.25)
    xlim([-1 3])
    set(gca,'XTick',[],'YTick',[],'XColor', 'none','YColor','none')
    if iSVD==1
        legend(handles,params.triallabels(1:2),'Location','NorthEast'); legend boxoff;
        sc = scalebar;
        sc.XLen = 1;
        sc.XUnit = 'sec';
        sc.YLen = 0;
        sc.YUnit = '';
        sc.Position = [-0.5, curylim(2)/2];  %move the whole SCALE position.
    end
end

%Audio: 
for iSVD = 1:params.nSVDs
    subplot(1,6,iSVD+3); hold all;
    handles = [];
    for iSplit = 1:2
%         handles(iSplit) = plot(params.xtime_video*1e-6,squeeze(meantoplot(:,iSVD,iSplit+2)),'Color',params.colors_splits{iSplit+2},'LineWidth',1);
        
        h = shadedErrorBar(params.xtime_video*1e-6,squeeze(meantoplot(:,iSVD,iSplit+2)),squeeze(errortoplot(:,iSVD,iSplit+2)),{'-','markerfacecolor',params.colors_splits{iSplit},'LineWidth',1},params.transparent);
        h.mainLine.Color = params.colors_splits{iSplit};    h.patch.FaceColor = params.colors_splits{iSplit};
        delete(h.edge(1:2));
        handles(iSplit) = h.mainLine;
    end
    %Figure make up:
    if iSVD==2
        title('Auditory')
    end
    curylim = get(gca,'ylim');
    plot([0 0],curylim,'-','Color',[0.6 0.6 0.6],'LineWidth',0.25)
    xlim([-1 3])
    set(gca,'XTick',[],'YTick',[],'XColor', 'none','YColor','none')
    if iSVD==1
        legend(handles,params.triallabels(3:4),'Location','NorthEast'); legend boxoff;
    end
end

filename = sprintf('motSVD_diffStim_%s.eps',params.exSes{1});
export_fig(fullfile(params.savedir,filename),gcf);


