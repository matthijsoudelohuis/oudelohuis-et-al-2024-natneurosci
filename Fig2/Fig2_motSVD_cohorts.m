%% Oude Lohuis et al. 2024 Nat Neurosci
% Triple dissociation of auditory, visual, and motor processing in primary visual cortex
% MOL (C) 2023

% Script analyzes motion in the video for the three cohorts

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

params.smoothSVD            = 0;
params.nSVDs                = 34;

params.savedir              = 'E:\Documents\PhD\Figures\Project CHDET\Results - auV1\13MovieBehavior';

%% Get Data:
[Data] = MOL_GetData('E:','CHDET',params.Experiments,{},{},{'sessionData' 'trialData_newtrials' 'spikeData' 'videoData'});
sessionData     = Data.sessionData;
trialData       = Data.trialData;
spikeData       = Data.spikeData;
videoData       = Data.videoData;

%% Remove last 20 trials:
trialData = MOL_RemoveLastnTrials(trialData,20);

%% Remove sessions that are passive/active from 2044 and 2045s: 
sesids                      = sessionData.session_ID(~ismember(sessionData.mousename,{'2044' '2045'}));
fprintf('Removed %d/%d sessions from 2044 and 2045 with active and passive epochs\n',length(sessionData.session_ID)-length(sesids),length(sessionData.session_ID));
[sessionData, trialData, spikeData, videoData]        = MOL_getTempPerSes(sesids,sessionData,trialData,spikeData, videoData);

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
idx             = ismember(spikeData.area,params.area);
fprintf('Filtered %d/%d neurons based on area\n',sum(idx),length(spikeData.session_ID));
for iF = 1:length(spikeFields)
    spikeData.(spikeFields{iF}) = spikeData.(spikeFields{iF})(idx);
end

%% Only sessions with video:
% Filter sessions based on containing motSVD data:
sesids = videoData.session_ID(~cellfun(@isempty,videoData.motSVD));
fprintf('Removed %d/%d sessions without videoData\n',length(sessionData.session_ID)-length(sesids),length(sessionData.session_ID));
[sessionData,trialData,spikeData,videoData] = MOL_getTempPerSes(sesids,sessionData,trialData,spikeData,videoData);

%% Report dataset:
fprintf('Dataset: %d sessions, %d trials, %d neurons, %d videos\n',length(sessionData.session_ID),length(trialData.session_ID),length(spikeData.session_ID),length(videoData.session_ID));
nSessions = length(sessionData.session_ID);

%% Save dataset:
save('Dataset2_1.mat','params','sessionData','trialData','spikeData','videoData','-v7.3')

%% Or load dataset
load Dataset2_1.mat





%% Make psth of mot svd from video: (takes some time)
[hist_mat,hist_mat_tot,hist_mat_z,params] = MOL_PSTH_video(sessionData,trialData,videoData,params);

nTrials             = length(trialData.stimChange);

%% Fig S3e: mean motion response to the different trial types for an example MST session:

[tempsessionData,temptrialData,tempspikeData,tempvideoData] = MOL_getTempPerSes({'20302020012158'},sessionData,trialData,spikeData,videoData);
nSVDs       = 3; %amount of pca's to plot

if params.smoothSVD
    for iSVD = 1:nSVDs
        tempvideoData.motSVD{1}(:,iSVD) = smooth(tempvideoData.motSVD{1}(:,iSVD),params.smoothSVD,'sgolay');
    end
end

nTempTrials         = length(temptrialData.stimChange);

hist_mat_exses      = NaN(params.nTimebins_video,nSVDs,nTempTrials);

fprintf('Computing svd response for trial        \n');
for iTrial = 1:nTempTrials
    fprintf(repmat('\b', 1, numel([num2str(iTrial-1) num2str(nTempTrials)])+1));
    fprintf('%d/%d',iTrial,nTempTrials);
    idx                                         = tempvideoData.ts{1}>temptrialData.(params.AlignOn)(iTrial)+params.t_pre & tempvideoData.ts{1}<temptrialData.(params.AlignOn)(iTrial)+params.t_post;
    hist_mat_exses(1:sum(idx),:,iTrial)     	= tempvideoData.(params.videofield){1}(idx,1:nSVDs);
end
fprintf('\n')

splits          = {};
splits{1}       = strcmp(temptrialData.trialType,'X') & temptrialData.vecResponse==3;
splits{2}       = strcmp(temptrialData.trialType,'X') & temptrialData.vecResponse==2;
splits{3}       = strcmp(temptrialData.trialType,'Y') & temptrialData.vecResponse==3;
splits{4}       = strcmp(temptrialData.trialType,'Y') & temptrialData.vecResponse==1;
nSplits         = length(splits);

params.triallabels  = {'Visual Miss' 'Visual Hit' 'Auditory Miss' 'Auditory Hit'};
params.trialcolors  = {[0.25 0.25 1] [0 0 0.8] [1 0.5 0.5] [0.8 0 0]};
params.triallines   = {'-.' '-' '-.' '-'};

examptoplot = NaN(length(splits),nSVDs,params.nTimebins_video);

%Select random trials from categories:
for iSplit = 1:nSplits
    idx = find(splits{iSplit});
    params.exampleTrials(iSplit) = temptrialData.trial_ID(idx(randi(numel(idx))));
end
%Or select these example trials:
params.exampleTrials = {'20302020012158350'   '20302020012158399'   '20302020012158289'    '20302020012158032'};

for iSVD = 1:nSVDs
    for iSplit = 1:nSplits
        %         examptoplot(iSplit,iSVD,:)  = hist_mat(:,iSVD,strcmp(trialData.trial_ID,params.exampleTrials{iSplit}));
        examptoplot(iSplit,iSVD,:)  = hist_mat_exses(:,iSVD,strcmp(temptrialData.trial_ID,params.exampleTrials{iSplit}));
    end
end

meantoplot = NaN(nSplits,nSVDs,params.nTimebins_video);
errortoplot = NaN(nSplits,nSVDs,params.nTimebins_video);
for iSVD = 1:nSVDs
    for iSplit = 1:nSplits
        meantoplot(iSplit,iSVD,:)    = nanmean(hist_mat_exses(:,iSVD,splits{iSplit}),3);
        errortoplot(iSplit,iSVD,:)   = nanstd(hist_mat_exses(:,iSVD,splits{iSplit},:),[],3);
        errortoplot(iSplit,iSVD,:)   = errortoplot(iSplit,iSVD,:) / sqrt(sum(splits{iSplit}));
    end
end

handles = [];
figure; hold all; set(gcf,'units','normalized','Position',[0.15 0.36 0.3 0.27],'color','w');
for iSVD = 1:nSVDs
    for iSplit = 1:nSplits
        subplot(nSVDs,nSplits+1,(iSVD-1)*(nSplits+1)+iSplit); hold all;
        plot(params.xtime_video*1e-6,squeeze(examptoplot(iSplit,iSVD,:)),'Color',params.trialcolors{iSplit},'LineWidth',1)
        xlim(params.xtime_video([1 end])*1e-6)
        ylim([-1000 1000])
        plot([0 0],[-500 500],'-','Color',[0.6 0.6 0.6],'LineWidth',0.25)
        %         set(gca,'XColor', 'none','YColor','none')
        set(gca,'xtick',[],'ytick',[])
        %         axis off
        if iSplit==1
            ylabel(sprintf('PCA %d',iSVD))
        end
    end
end

for iSVD = 1:nSVDs
    subplot(nSVDs,nSplits+1,(iSVD-1)*(nSplits+1)+5); hold all;
    for iSplit = 1:length(splits)
        handles(iSplit) = plot(params.xtime_video*1e-6,squeeze(meantoplot(iSplit,iSVD,:)),'Color',params.trialcolors{iSplit},'LineWidth',1);
    end
    xlim(params.xtime_video([1 end-1])*1e-6)
    ylim([-1000 1000])
    set(gca,'XColor', 'none','YColor','none')
    plot([0 0],[-500 500],'-','Color',[0.6 0.6 0.6],'LineWidth',0.25)
end

xlim([params.xtime_video(1)*1e-6 params.xtime_video(end)*1e-6])
legend(handles,params.triallabels,'Location','SouthEast'); legend boxoff;

tightfig;

filename = sprintf('motSVD_exSession.eps');
% export_fig(fullfile(params.savedir,filename),gcf);

%% Fig 2b: mean motion response for the different cohorts:

params.nSplits  = 4;
splits          = {};
splits{1}       = strcmp(trialData.trialType,'X') & trialData.visualOriChangeNorm==2;
splits{2}       = strcmp(trialData.trialType,'X') & trialData.visualOriChangeNorm==3;
splits{3}       = strcmp(trialData.trialType,'Y') & trialData.audioFreqChangeNorm==2;
splits{4}       = strcmp(trialData.trialType,'Y') & trialData.audioFreqChangeNorm==3;

params.triallabels = {'Vthr' 'Vmax' 'Athr' 'Amax'};

meantoplot      = NaN(params.nExperiments,params.nSplits,params.nTimebins_video); %Init matrix vars
errortoplot     = NaN(params.nExperiments,params.nSplits,params.nTimebins_video);

for iExp = 1:params.nExperiments %Loop over cohorts
    idx_exp                     = ismember(trialData.session_ID,sessionData.session_ID(strcmp(sessionData.Experiment,params.Experiments(iExp)))); %Get all trials from specific animals (NE, UST, MST)
    for iSplit = 1:4
        idx_all                 = idx_exp & splits{iSplit};
        meantoplot(iExp,iSplit,:)      = nanmean(hist_mat_z(:,idx_all),2);
        errortoplot(iExp,iSplit,:)     = nanstd(hist_mat_z(:,idx_all),[],2) / sqrt(sum(idx_all));
    end
end

%Plot figure:
figure; hold all; set(gcf,'units','normalized','Position',[0.05 0.35 0.4 0.14],'color','w');
for iSplit = 1:4
    subplot(1,4,iSplit); hold all;
    handles = [];
    for iExp = 1:params.nExperiments
        handles(iExp) = plot(params.xtime_video*1e-3,squeeze(meantoplot(iExp,iSplit,:)),'Color',params.colors_experiments{iExp},'LineWidth',1);
    end
    
    %Figure make up:
    title(params.triallabels{iSplit})
    plot([0 0],[-0.5 1],'-','Color',[0.6 0.6 0.6],'LineWidth',0.25)
    xlim([-1e3 3e3])
    ylim([-0.2 3.5])
    set(gca,'XTick',[],'YTick',[],'XColor', 'none','YColor','none')
    if iSplit==1
        sc = scalebar;
        sc.XLen = 1000;
        sc.XUnit = 'sec';
        sc.YLen = 1;
        sc.YUnit = 'z-score';
        sc.Position = [-1, 1];  %move the whole SCALE position.
        legend(handles,params.ExperimentLabels,'Location','NorthEast'); legend boxoff;
    end
    if iSplit == 4
        plot([-0.2e3 0.3e3 0.3e3 -0.2e3 -0.2e3],[-0.2 -0.2 1 1 -0.2],'k:','LineWidth',0.5)
    end
    
end

filename = sprintf('motSVD_avgMot_3cohorts_trialtypes.eps');
export_fig(fullfile(params.savedir,filename),gcf);

%% Fig 2b: plot close up figure:
figure; hold all; set(gcf,'units','normalized','Position',[0.3 0.35 0.1 0.17],'color','w');
iSplit = 4;
handles = [];
for iExp = 1:params.nExperiments
    handles(iExp) = plot(params.xtime_video*1e-3,squeeze(meantoplot(iExp,iSplit,:)),'Color',params.colors_experiments{iExp},'LineWidth',1);
end

%Figure make up:
title(params.triallabels{iSplit})
plot([0 0],[-0.5 1],'-','Color',[0.6 0.6 0.6],'LineWidth',0.25)
xlim([-0.2e3 0.3e3])
ylim([-0.2 1])
set(gca,'XTick',(-1:0.05:2)*1e3,'YTick',[0 0.5 1 1.5],'XTickLabelRotation',45)

filename = sprintf('motSVD_avgMot_3cohorts_Amax.eps');
% export_fig(fullfile(params.savedir,filename),gcf);

%%











%% Neural correlates to video ME:



%% Fig. S3a: rate and video ME in example cells: 

params.cell_IDs = {'20092018082331259' '20292019121311297'};%One MST, one UST neuron %'20122018081611182' '20112018080911216'
selec = 5000:6000;
close all;

for iN = 1:length(params.cell_IDs)
    
    idx_cell        = strcmp(spikeData.cell_ID,params.cell_IDs{iN});
    sesid           = spikeData.session_ID(idx_cell);
    
    [tempsessionData,temptrialData,tempspikeData,tempvideoData] = MOL_getTempPerSes(sesid,sessionData,trialData,spikeData,videoData);

    xtime = tempvideoData.ts{1}(selec);
    
    spikerate   = histcounts(spikeData.ts{idx_cell},[xtime-20e3; xtime(end)+20e3]);
    video_ME    = sum(abs(tempvideoData.motSVD{1}(selec,:)),2);
    
    spikerate   = smooth(spikerate,7,'lowess');
    video_ME    = smooth(video_ME,5,'lowess');
    
    figure;  hold all; set(gcf,'units','normalized','Position',[0.05 0.4 0.83 0.07],'color','w');
    title(params.cell_IDs{iN})
    plot(xtime,spikerate / max(spikerate),'k-');
    plot(xtime,video_ME / max(video_ME),'r-');
    
    %Figure make up:
    xlim([xtime(1) xtime(end)])
    ylim([0 1])
    set(gca,'XTick',[],'YTick',[0 1],'XColor', 'none')
    sc = scalebar;
    sc.XLen = 1e6;
    sc.XUnit = 'sec';
    sc.YLen = 1;
    sc.YUnit = 'z-score';
    sc.Position = [xtime(100), 0.5];  %move the whole SCALE position.
    
    [r] = corr(spikerate,video_ME);
    fprintf('%s: r=%1.2f\n',params.cell_IDs{iN},r)
    
    filename = sprintf('Corr_trace_%s.eps',params.cell_IDs{iN});
%     export_fig(fullfile(params.savedir,filename),gcf);
end

%% Report quantification correlation: 
nNeurons                    = length(spikeData.ts);

rmat = NaN(nNeurons,1);
pmat = NaN(nNeurons,1);

for iNeuron = 1:nNeurons %Loop over all neurons:
    
    sesid           = spikeData.session_ID(iNeuron);
    
    [tempsessionData,temptrialData,tempspikeData,tempvideoData] = MOL_getTempPerSes(sesid,sessionData,trialData,spikeData,videoData);

    xtime = tempvideoData.ts{1};
    
    fprintf(repmat('\b', 1, numel([num2str(iNeuron-1) num2str(nNeurons)])+2));
    fprintf('%d/%d\n',iNeuron,nNeurons);
    
    spikerate   = histcounts(spikeData.ts{iNeuron},[xtime-20e3; xtime(end)+20e3]);
    video_ME    = sum(abs(tempvideoData.motSVD{1}(:,:)),2);
    
    [rmat(iNeuron),pmat(iNeuron)] = corr(spikerate',video_ME,'rows','complete');
    rtemp = NaN(1,1000);
    for i = 1:1000
        rtemp(i) = corr(spikerate',circshift(video_ME,randi(length(video_ME))),'rows','complete');
    end
    pmat(iNeuron) = sum(rmat(iNeuron)>abs(rtemp));
end

%% Report number of significantly correlated neurons
% Results may vary because of random circular shift in permutation test above
p = 0.01 / nNeurons; %bonferroni corr
%how many of circular shuffles have more negative or more positive actual correlation than shuffles
sigmat = pmat<(p*1000) | pmat>((1-p)*1000); 
fprintf('%2.1f%% (%2.0f/%2.0f neurons) had a significant correlation (P<0.05, circular shuffle test), with a mean absolute correlation of %1.2f\n',(sum(sigmat) / nNeurons)*100,sum(sigmat),nNeurons,mean(abs(rmat)))

%% Compute correlation of firing rate to video over time: 

params.nSplits  = 2;
splits          = {};
splits{1}       = strcmp(trialData.trialType,'X');
splits{2}       = strcmp(trialData.trialType,'Y');
params.labels_splits        = {'V' 'Vmiss' 'A' 'Amiss'};
params.labels_splits        = {'V' 'A'};

params.zscore               = 0;
params.smoothing            = 1;
params.minTrialCond         = 10; %Combined trial number minimum

fprintf('Computing average firing rate response for neuron        \n');

params.binsize              = 10e3;
%ReConstruct bin edges and time axis
params.edges                = [params.t_pre:params.binsize:params.t_post] - params.binsize/2;                    %#ok<NBRAK> Define edges of the bins
params.xtime                = params.t_pre:params.binsize:params.t_post-params.binsize;                        %Define time axis
params.nTimebins            = length(params.xtime); %number of time bins

nNeurons                    = length(spikeData.ts);
corrmat                     = NaN(nNeurons,params.nSplits,params.nTimebins);

for iNeuron = 1:nNeurons %Loop over all neurons:
    fprintf(repmat('\b', 1, numel([num2str(iNeuron-1) num2str(nNeurons)])+2));
    fprintf('%d/%d\n',iNeuron,nNeurons);
    
    sesid       = spikeData.session_ID(iNeuron);
    idx_ses     = strcmp(trialData.session_ID,sesid);
    %Compute histogram:
    events_ts                   = trialData.(params.AlignOn)(idx_ses);
    rate_mat                    = calc_psth(events_ts,spikeData.ts{iNeuron},params);    %Construct histogram matrix:
    
    for iSplit = 1:params.nSplits %Store the mean response for each of these splits
        idx_trial   = splits{iSplit}(idx_ses);
        
        temp_mot    = hist_mat_tot(:,idx_ses & splits{iSplit});
        temp_rate   = rate_mat(idx_trial,:)';
        
        temp_mot   = interp1(params.xtime_video,temp_mot,params.xtime);
        
        for iBin = 1:params.nTimebins
            if sum(temp_rate(iBin,:)~=0)>params.minTrialCond
                corrmat(iNeuron,iSplit,iBin) = corr(temp_mot(iBin,:)',temp_rate(iBin,:)','rows','complete');
            end
        end
    end
end

%% Fig 2f: correlation over auditory or visual trials:
meantoplot                  = squeeze(nanmean(corrmat,1));
errortoplot                 = squeeze(nanstd(corrmat,[],1)) / sqrt(nNeurons);

params.trialcolors = {[0 0 0.8000] [0.8000 0 0]};
figure;  hold all; set(gcf,'units','normalized','Position',[0.05 0.4 0.18 0.24],'color','w');

for iSplit = [1 2] %loop across modalities
    h = shadedErrorBar(params.xtime,squeeze(meantoplot(iSplit,:)),squeeze(errortoplot(iSplit,:)),{'-','Color',params.trialcolors{iSplit},'LineWidth',1},0);
    delete(h.edge(1:2));
end
xlim([-0.2e6 0.6e6])
ylim([0.02 0.12])
set(gca,'XTick',-4e6:0.1e6:4e6,'XTickLabels',(-4:0.1:4)*1e3)
grid on
filename = sprintf('Corr_Temporal_MotSVD_3cohorts_AVtrialtypes.eps');
% export_fig(fullfile(params.savedir,filename),gcf); %Doesnt work with shaded error, save manually

