%% Script that analyzes movement across active and passive task versions:
% MOL (C) 2021
startover

%% 
params                      = struct();
params                      = params_histresponse_auV1(params);
params.AlignOn              = 'stimChange';         %On which timestamp to align trials to

params.videofield           = 'motSVD';

params.fs                   = 25; %Hz

params.t_pre                = -1e6;
params.t_post               = 4e6;

params                      = MOL_getColors_CHDET(params);

params.area                 = 'V1';

params.savedir              = 'E:\Documents\PhD\Figures\Project CHDET\Results - auV1\3ActivePassive';

%% Get data:
[Data] = MOL_GetData('E:','CHDET','ChangeDetectionConflict',{'2044' '2045'},[],{'sessionData' 'trialData' 'spikeData' 'videoData'});
sessionData     = Data.sessionData;
trialData       = Data.trialData;
spikeData       = Data.spikeData;
videoData       = Data.videoData;

%% Remove last 1 trials:
trialData           = MOL_RemoveLastnTrials(trialData,1);

%% Get sessions with Engaged and passive part:
sesids = {    '20442021042434'    '20442021042435'    '20442021042837'    '20442021042838'    '20452021042849'    '20452021042850'...
    '20452021042952'    '20452021042953'    '20452021043055'    '20452021043056' '20312020012205'};
[sessionData,trialData,spikeData,videoData] = MOL_getTempPerSes(sesids(:),sessionData,trialData,spikeData,videoData);

%% Get sessions with video:
sesids = videoData.session_ID(~cellfun(@isempty,videoData.motSVD));
[sessionData,trialData,spikeData,videoData] = MOL_getTempPerSes(sesids(:),sessionData,trialData,spikeData,videoData);

%% Set trialData to active and passive:
nSessions       = length(sessionData.session_ID);
nTrials         = length(trialData.session_ID);
trialData.State = NaN(nTrials,1);

for iSes = 1:nSessions
    idx = strcmp(trialData.session_ID,sessionData.session_ID(iSes));
    if strcmp(sessionData.State(iSes),'Awake')
        trialData.State(idx)        = 1;
    elseif strcmp(sessionData.State(iSes),'Behaving')
        trialData.State(idx)        = 2;
    end
end

%% Manually set periods of removed lick spout in session of 31 to passive:
idx = strcmp(trialData.session_ID,'20312020012205') & ismember(trialData.trialNum,[53:104 164:215 268:325]);
% idx = strcmp(trialData.session_ID,'20442021042128') & ismember(trialData.trialNum,[53:104 164:215 268:325]);
trialData.State(idx)         = 1;

%% Combine trialData from active and passive part of same session:
sesids = {    '20442021042434'    '20442021042435';
    '20442021042837'    '20442021042838';
    '20452021042849'    '20452021042850';
    '20452021042952'    '20452021042953';
    '20452021043055'    '20452021043056'};

nTrials         = length(trialData.session_ID);
for iTrial = 1:nTrials
    temp = strcmp(trialData.session_ID(iTrial),sesids(:,2));
    if any(temp)
        trialData.session_ID(iTrial) = sesids(temp,1);
    end
end

%% Combine neurons:
nNeurons = length(spikeData.ts);
for iNeuron = 1:nNeurons %Loop over all neurons:
    temp = strcmp(spikeData.session_ID(iNeuron),sesids(:,2));
    if any(temp)
        spikeData.session_ID(iNeuron) = sesids(temp,1);
    end
end
spikeFields     = fieldnames(spikeData);
[C,IA,~] = unique(spikeData.cell_ID);

for iN = 1:length(C)
    dupl = find(strcmp(spikeData.cell_ID,C(iN)) & ~ismember(1:nNeurons,IA)');
    if dupl
        spikeData.ts{IA(iN)} = [spikeData.ts{IA(iN)} spikeData.ts{dupl}];
    end
end

for iF = 1:length(spikeFields)
    spikeData.(spikeFields{iF}) = spikeData.(spikeFields{iF})(IA);
end
fprintf('Merged spikeData\n');

%% Combine videos:
videoFields =  fieldnames(videoData);
videoFields = videoFields(~ismember(videoFields,{'fs','session_ID'}));

for iSes = 1:size(sesids,1)
    idx_1 = strcmp(sessionData.session_ID,sesids{iSes,1});
    if any(idx_1)
        idx_2 = strcmp(sessionData.session_ID,sesids{iSes,2});
        for iF = 1:length(videoFields)
            videoData.(videoFields{iF}){idx_1} = [videoData.(videoFields{iF}){idx_1}; videoData.(videoFields{iF}){idx_2}];
        end
    end
end

fprintf('Merged videoData\n');

%% Filter out duplicate session and video data:
sesids = sessionData.session_ID(~ismember(sessionData.session_ID,sesids(:,2)));
[sessionData,trialData,spikeData,videoData] = MOL_getTempPerSes(sesids,sessionData,trialData,spikeData,videoData);
nSessions       = length(sessionData.session_ID);
nTrials         = length(trialData.session_ID);

%% Make psth of video mot SVD relative to stims:
params.t_pre        = -1e6;
params.t_post       = 4e6;
params.nSVDs        = 25;

[video_hist,video_hist_tot,video_hist_z,params] = MOL_PSTH_video(sessionData,trialData,videoData,params);

%% Make figure of mean motion response for the different states:

params.nSplits  = 2;
splits          = {};
splits{1}       = strcmp(trialData.trialType,'Y') & trialData.State==1;
splits{2}       = strcmp(trialData.trialType,'Y') & trialData.State==2;

meantoplot      = NaN(params.nSplits,params.nTimebins_video); %Init matrix vars
errortoplot     = NaN(params.nSplits,params.nTimebins_video);

for iSplit = 1:params.nSplits 
    idx_all                         = splits{iSplit};
    meantoplot(iSplit,:)            = nanmean(video_hist_z(:,idx_all),2);
    errortoplot(iSplit,:)           = nanstd(video_hist_z(:,idx_all),[],2) / sqrt(sum(idx_all));
end

%Plot figure:
figure; hold all; set(gcf,'units','normalized','Position',[0.3 0.45 0.23 0.3],'color','w');
    handles = [];
for iSplit = 1:params.nSplits
    handles(iSplit) = plot(params.xtime_video*1e-6,squeeze(meantoplot(iSplit,:)),'Color',params.colors_state{iSplit},'LineWidth',1);
end

%Figure make up:
% title(params.triallabels{iSplit})
plot([0 0],[-0.5 1],'-','Color',[0.6 0.6 0.6],'LineWidth',0.25)
xlim([-0.5 2])
ylim([-0.2 1])
% set(gca,'XTick',[],'YTick',[],'XColor', 'none','YColor','none')
% sc = scalebar;
% sc.XLen = 0.5;
% sc.XUnit = 'sec';
% sc.YLen = 0.5;
% sc.YUnit = 'z-score';
% sc.Position = [-0.5, 0.5];  %move the whole SCALE position.
legend(handles,params.labels_state,'Location','NorthEast'); legend boxoff;

filename = sprintf('motSVD_ActPas_Autrials.eps');
export_fig(fullfile(params.savedir,filename),gcf);

%% Make figure of mean motion response for the different states:

params.nSplits  = 2;
splits          = {};
splits{1}       = strcmp(trialData.trialType,'X') & trialData.State==1;
splits{2}       = strcmp(trialData.trialType,'X') & trialData.State==2;

meantoplot      = NaN(params.nSplits,params.nTimebins_video); %Init matrix vars
errortoplot     = NaN(params.nSplits,params.nTimebins_video);

for iSplit = 1:params.nSplits 
    idx_all                         = splits{iSplit};
    meantoplot(iSplit,:)            = nanmean(video_hist_z(:,idx_all),2);
    errortoplot(iSplit,:)           = nanstd(video_hist_z(:,idx_all),[],2) / sqrt(sum(idx_all));
end

%Plot figure:
figure; hold all; set(gcf,'units','normalized','Position',[0.3 0.45 0.23 0.3],'color','w');
    handles = [];
for iSplit = 1:params.nSplits
    handles(iSplit) = plot(params.xtime_video*1e-6,squeeze(meantoplot(iSplit,:)),'Color',params.colors_state{iSplit},'LineWidth',1);
end

%Figure make up:
% title(params.triallabels{iSplit})
plot([0 0],[-0.5 1],'-','Color',[0.6 0.6 0.6],'LineWidth',0.25)
xlim([-0.5 2])
ylim([-0.2 1])
% set(gca,'XTick',[],'YTick',[],'XColor', 'none','YColor','none')
% sc = scalebar;
% sc.XLen = 0.5;
% sc.XUnit = 'sec';
% sc.YLen = 0.5;
% sc.YUnit = 'z-score';
% sc.Position = [-0.5, 0.5];  %move the whole SCALE position.
legend(handles,params.labels_state,'Location','NorthEast'); legend boxoff;

filename = sprintf('motSVD_ActPas_Vistrials.eps');
export_fig(fullfile(params.savedir,filename),gcf);

%% Make figure of mean motion response for the different states:

params.nSplits  = 2;
params.triallabels = {'Pas' 'Act'};

meantoplot      = NaN(params.nSplits,params.nTimebins_video); %Init matrix vars
errortoplot     = NaN(params.nSplits,params.nTimebins_video);

for iSplit = 1:params.nSplits 
    idx_all                         = splits{iSplit};
    meantoplot(iSplit,:)            = nanmean(video_hist_z(:,idx_all),2);
    errortoplot(iSplit,:)           = nanstd(video_hist_z(:,idx_all),[],2) / sqrt(sum(idx_all));
end

%Plot figure:
figure; hold all; set(gcf,'units','normalized','Position',[0.3 0.45 0.23 0.3],'color','w');
handles = [];

for iSes = 1:nSessions
    idx                 = strcmp(trialData.session_ID,sessionData.session_ID(iSes)) & strcmp(trialData.trialType,'Y');
    meantoplot          = nanmean(video_hist_z(:,idx),2);
%     iSplit = sessionData.State{iSes}-1;
    iSplit = strcmp(sessionData.State{iSes},'Behaving')+1;
%     handles(iSes) = plot(params.xtime_video*1e-6,meantoplot,'Color',params.colors_state{iSplit},'LineWidth',1);
    handles = plot(params.xtime_video*1e-6,meantoplot,'Color',params.colors_state{iSplit},'LineWidth',1);
    
end

% 
% for iSplit = 1:params.nSplits
%     handles(iSplit) = plot(params.xtime_video*1e-6,squeeze(meantoplot(iSplit,:)),'Color',params.colors_state{iSplit},'LineWidth',1);
% end

%Figure make up:
% title(params.triallabels{iSplit})
plot([0 0],[-0.5 1],'-','Color',[0.6 0.6 0.6],'LineWidth',0.25)
xlim([-0.5 2])
ylim([-0.2 1])
set(gca,'XTick',[],'YTick',[],'XColor', 'none','YColor','none')
sc = scalebar;
sc.XLen = 0.5;
sc.XUnit = 'sec';
sc.YLen = 0.5;
sc.YUnit = 'z-score';
sc.Position = [-0.5, 0.5];  %move the whole SCALE position.
legend(handles,params.labels_state,'Location','NorthEast'); legend boxoff;

% filename = sprintf('motSVD_avgMot_3cohorts_trialtypes.eps');
% export_fig(fullfile(params.savedir,filename),gcf);
