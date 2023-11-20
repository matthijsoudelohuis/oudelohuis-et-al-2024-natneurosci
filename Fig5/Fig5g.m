%% Oude Lohuis et al. 2023 Nat Neurosci
% Triple dissociation of auditory, visual, and motor processing in primary visual cortex
% MOL (C) 2023

% Script analyzes video motion in sessions with and without muscimol

startover

%% Parameter settings:
params                      = params_histresponse_auV1;
params.AlignOn              = 'stimChange';         %On which timestamp to align trials to

params.videofield           = 'motSVD';

params.area                 = 'V1';

params.fs                   = 25; %Hz

params.t_pre                = -1e6;
params.t_post               = 4e6;

params.smoothSVD            = 0;
params.nSVDs                = 25;
params.Project              = 'CHDET'; %Which versions of the task to load data from
params.Experiment           = 'ChangeDetectionConflict'; %Which versions of the task to load data from
params.savedir              = 'E:\Documents\PhD\Figures\Project CHDET\Results - auV1\7MuscimolA1';

params                      = MOL_getColors_CHDET(params);

%% Configurations from excel overview of recordings per project:
RootDir                 = 'E:\';
MainDataDir             = fullfile(RootDir,'Data',params.Project);
[~,~,RawLibrary]        = xlsread(fullfile(RootDir,'Data',params.Project,sprintf('%s_Sessions_Overview.xlsx',params.Project)));
RawLibrary              = RawLibrary(~cellfun(@isnumeric,RawLibrary(:,1)),:); %Trim the RawLibrary to entered values
[XlsRows,XlsColumns]    = size(RawLibrary);

ExpIdx                  = strcmp(RawLibrary(:,strcmp(RawLibrary(1,:),'MuscimolArea')),'A1');
Sesselec                = RawLibrary(ExpIdx,strcmp(RawLibrary(1,:),'Rec_datetime'));

ExpIdx                  = strcmp(RawLibrary(:,strcmp(RawLibrary(1,:),'MuscimolArea')),'A1') | strcmp(RawLibrary(:,strcmp(RawLibrary(1,:),'SalineArea')),'A1');
Sesselec                = RawLibrary(ExpIdx,strcmp(RawLibrary(1,:),'Rec_datetime'));

%% Get data:
[Data] = MOL_GetData('E:','CHDET',params.Experiment,{},Sesselec,{'sessionData' 'trialData_newtrials' 'spikeData' 'videoData'});
sessionData     = Data.sessionData;
trialData       = Data.trialData;
spikeData       = Data.spikeData;
videoData       = Data.videoData;

%% Remove last 20 trials:
trialData = MOL_RemoveLastnTrials(trialData,20);

%% Filter sessions based on containing motSVD data:
sesids = videoData.session_ID(~cellfun(@isempty,videoData.motSVD));
fprintf('Removed %d/%d sessions without videoData\n',length(sessionData.session_ID)-length(sesids),length(sessionData.session_ID));
[sessionData,trialData,spikeData,videoData] = MOL_getTempPerSes(sesids,sessionData,trialData,spikeData,videoData);

%% Save dataset:
save('Dataset5_4.mat','params','sessionData','trialData','spikeData','videoData')

%% Or start script from saved dataset:
load Dataset5_4.mat

%% Report dataset:
fprintf('Dataset: %d sessions, %d trials, %d neurons, %d videos\n',length(sessionData.session_ID),length(trialData.session_ID),length(spikeData.session_ID),length(videoData.session_ID));
nSessions = length(sessionData.session_ID);



%% Make psth of mot svd from video:
[hist_mat,hist_mat_tot,hist_mat_z,params] = MOL_PSTH_video(sessionData,trialData,videoData,params);

nTrials             = length(trialData.stimChange);

%% Fig 5g: mean motion response

params.lines_mans           = {'-' ':'};

params.nSplits  = 1;
splits          = {};
splits{1}       = strcmp(trialData.trialType,'Y') & trialData.vecResponse==1;

% params.triallabels = {'Vthr' 'Vmax' 'Athr' 'Amax'};
params.triallabels = {'Amax'};

meantoplot      = NaN(2,params.nSplits,params.nTimebins_video); %Init matrix vars
errortoplot     = NaN(2,params.nSplits,params.nTimebins_video);

idx_man = [];
idx_man(:,1)    = ismember(trialData.session_ID,sessionData.session_ID(strcmp(sessionData.SalineArea,'A1')));
idx_man(:,2)    = ismember(trialData.session_ID,sessionData.session_ID(strcmp(sessionData.MuscimolArea,'A1')));

for iMan = 1:2
    for iSplit = 1:params.nSplits
        idx_all                 = idx_man(:,iMan) & splits{iSplit};
        meantoplot(iMan,iSplit,:)      = nanmean(hist_mat_z(:,idx_all),2);
        errortoplot(iMan,iSplit,:)     = nanstd(hist_mat_z(:,idx_all),[],2) / sqrt(sum(idx_all));
    end
end

%Plot figure:
figure; hold all; set(gcf,'units','normalized','Position',[0.05 0.35 0.18 0.19],'color','w');
for iSplit = 1:params.nSplits
    subplot(1,params.nSplits,iSplit); hold all;
    handles = [];
    for iMan = 1:2
        %         handles(iExp) = plot(params.xtime_video*1e-3,squeeze(meantoplot(iSplit,:)),'Color',params.colors_experiments{iExp},'LineWidth',1);
        h = shadedErrorBar(params.xtime_video*1e-3,squeeze(meantoplot(iMan,iSplit,:)),squeeze(errortoplot(iMan,iSplit,:)),{'k','Color',params.colors_ztrials{5-iMan},'LineWidth',1,'LineStyle',params.lines_mans{iMan}},0);
        delete(h.edge(:));
        handles(iMan) = h.mainLine;
    end
    %Figure make up:
%     title(params.triallabels{iSplit})
    plot([0 0],[-0.5 1],'-','Color',[0.6 0.6 0.6],'LineWidth',0.25)
    xlim([-1e3 3e3])
    ylim([-0.6 4])
    set(gca,'XTick',[],'YTick',[],'XColor', 'none','YColor','none')
    if iSplit==1
        sc = scalebar;
        sc.XLen = 1000;
        sc.XUnit = 'sec';
        sc.YLen = 1;
        sc.YUnit = 'z-score';
        sc.Position = [-1, 1];  %move the whole SCALE position.
        legend(handles,{'Sal' 'Mus'},'Location','NorthEast'); legend boxoff;
    end
end

filename = sprintf('motSVD_avgMot_musc.eps');
export_fig(fullfile(params.savedir,filename),gcf);

