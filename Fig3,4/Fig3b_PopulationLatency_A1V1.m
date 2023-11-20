%% Oude Lohuis et al. 2024 Nat Neurosci
% Triple dissociation of auditory, visual, and motor processing in primary visual cortex
% MOL (C) 2023

startover

%% Parameter settings:
params                      = params_histresponse(); % All time is in microseconds
params.Experiments          = {'ChangeDetectionConflictDecor' 'VisOnlyTwolevels' 'ChangeDetectionConflict' }; %Which versions of the task to load data from
params.ExperimentLabels     = {'NE' 'UST' 'MST'}; %Which versions of the task to load data from

params.AlignOn              = 'stimChange';         %On which timestamp to align trials to

params.minNneurons          = 10;
params.zthreshold           = 1;

params.areas                = {'A1' 'V1'};
params.nAreas               = length(params.areas);

params.savedir              = 'E:\Documents\PhD\Figures\Project CHDET\Results - auV1\19AC\';

%% Get data:
[Data] = MOL_GetData('E:','CHDET',params.Experiments,{},[],{'sessionData' 'trialData_newtrials' 'spikeData'});
sessionData     = Data.sessionData;
trialData       = Data.trialData;
spikeData       = Data.spikeData;

%% Remove last 20 trials:
trialData = MOL_RemoveLastnTrials(trialData,20);

%% Filter out neurons based on quality:
spikeData       = MOL_filterNeurons(sessionData,trialData,spikeData);

%% Filter out sessions with muscimol:
sesids              = sessionData.session_ID(cellfun(@isempty,sessionData.MuscimolArea));
fprintf('Removed %d/%d sessions  with muscimol manipulations\n',length(sessionData.session_ID)-length(sesids),length(sessionData.session_ID));
[sessionData, trialData, spikeData]        = MOL_getTempPerSes(sesids,sessionData,trialData,spikeData);

%% Filter out trials with photostimulation in V1:
sesids              = sessionData.session_ID(sessionData.UseOpto & strcmp(sessionData.PhotostimArea,'V1'));
idx                 = trialData.hasphotostim==1 & ismember(trialData.session_ID,sesids);
trialFields         = fieldnames(trialData);
fprintf('Removed %d/%d trials  with photostim in V1\n',sum(idx),length(trialData.session_ID));
for iF = 1:length(trialFields)
    trialData.(trialFields{iF}) = trialData.(trialFields{iF})(~idx);
end

%% Set recordings in animal 1012 to be from LM based on histology:
idx = strcmp(spikeData.session_ID,'10122019041126') & strcmp(spikeData.area,'A1');
spikeData.area(idx) = deal({'LM'});

%% Filter neurons in area:
idx             = ismember(spikeData.area,params.areas);
spikeFields     = fieldnames(spikeData);
for iF = 1:length(spikeFields)
    spikeData.(spikeFields{iF}) = spikeData.(spikeFields{iF})(idx);
end

%% Get only sessions with neurons recorded:
[sessionData,trialData,spikeData] = MOL_getTempPerSes(unique(spikeData.session_ID),sessionData,trialData,spikeData);
fprintf('Dataset: %d sessions, %d trials, %d neurons\n',length(sessionData.session_ID),length(trialData.session_ID),length(spikeData.session_ID));

%% Save dataset:
save('Dataset3_1.mat','params','sessionData','trialData','spikeData')

%% Or start script from saved dataset:
load Dataset3_1.mat




%% histogram parameters:
params.zscore               = 0;
params.smoothing            = 1;
params.binsize              = 1e3;             %Size of the bins
params.conv_win             = 'chg';       %Type of window used for smoothing {flat, gaussian)
params.conv_twin            = 0.3e6;            %Window size for smoothing
params.conv_sigma           = 10e3;           %sd of gaussian window for smoothing

%Construct bin edges and time axis
params.edges                = [params.t_pre:params.binsize:params.t_post] - params.binsize/2;                    %#ok<NBRAK> Define edges of the bins
params.xtime                = params.t_pre:params.binsize:params.t_post-params.binsize;                        %Define time axis
params.nTimebins            = length(params.xtime); %number of time bins

%% For each neuron compute firing rate:s
fprintf('Computing responses for session     \n');

params.nSplits          = 2;
nNeurons                = length(spikeData.ts);
lastsesid               = []; %memory var to keep track of whether neuron comes from different session and design matrix needs to be reconstructed
fprintf('Computing average firing rate for neuron        \n');
zmat                    = NaN(nNeurons,params.nTimebins,params.nSplits);
params.minTrialCond     = 10;

for iNeuron = 1:nNeurons %Loop over all neurons:
    fprintf(repmat('\b', 1, numel([num2str(iNeuron-1) num2str(nNeurons)])+2));
    fprintf('%d/%d\n',iNeuron,nNeurons);
    
    if ~strcmp(lastsesid,spikeData.session_ID(iNeuron)) %construct new predictor matrix if neuron comes from a new session:
        %Get the relevant data for each session individually:
        temptrialData        = MOL_getTempPerSes(spikeData.session_ID(iNeuron),trialData);
        lastsesid            = spikeData.session_ID(iNeuron); %save this session_ID
    end
    
    %Compute histogram:
    events_ts               = temptrialData.(params.AlignOn);
    hist_mat                = calc_psth(events_ts,spikeData.ts{iNeuron},params);    %Construct histogram matrix:
    
    splits                  = {};
    splits{1}               = ismember(temptrialData.trialType,{'X'}) & temptrialData.visualOriChangeNorm==3;
    splits{2}               = ismember(temptrialData.trialType,{'Y'}) & temptrialData.audioFreqChangeNorm==3;
    
    for iSplit = 1:params.nSplits %Store the mean response for each of these splits
        if sum(splits{iSplit})>=params.minTrialCond
            zmat(iNeuron,:,iSplit) = nanmean(hist_mat(splits{iSplit},:),1);
        end
    end
end

%% For each session and each area take the mean firing rate across all neurons:
nSessions                   = length(sessionData.session_ID);
sesmat                      = NaN(nSessions,params.nTimebins,params.nSplits,params.nAreas);

for iSes = 1:nSessions %Loop over all neurons:
    for iArea = 1:2
        idx = strcmp(spikeData.session_ID,sessionData.session_ID(iSes)) & strcmp(spikeData.area,params.areas(iArea));
        
        if sum(idx)>=params.minNneurons
            sesmat(iSes,:,:,iArea) = nanmean(zmat(idx,:,:),1);
        end
    end
end

%% Zscore the population firing rate:
idx_bsl = params.xtime<-0.1e6;

for iSes = 1:nSessions %Loop over all neurons:
    for iSplit = 1:params.nSplits %Store the mean response for each of these splits
        for iArea = 1:2
%             sesmat(iSes,:,iSplit,iArea) = zscore(sesmat(iSes,:,iSplit,iArea));
            temp                        = sesmat(iSes,:,iSplit,iArea);
            sesmat(iSes,:,iSplit,iArea) = (temp - mean(temp(idx_bsl))) / std(temp(idx_bsl));
        end
    end
end

%% Compute latency to cross z-scored threshold, averaged across sessions:
poplat = NaN(2,2); %2x2 areas by modalities
for iArea = 1:2
    for iMod = 1:2
        poplat(iArea,iMod) = params.xtime(find(squeeze(nanmean(sesmat(:,:,iMod,iArea),1))>params.zthreshold,1))*1e-3;
    end
end

%% Extra figure to show individual sessions: (not in MS)
figure; hold all; set(gcf,'units','normalized','Position',[0.12 0.5 0.5 0.32],'color','w');
subplot(1,2,1); hold all;

for iSes = 1:nSessions %Loop over all sessions:
    plot(params.xtime,sesmat(iSes,:,1,1),'r','LineWidth',0.25)
    plot(params.xtime,sesmat(iSes,:,1,2),'b','LineWidth',0.25)
end
%Figure make up:
title('Visual Stim')
legend(params.areas)
xlim([-0.5e6 1e6])

subplot(1,2,2); hold all;
for iSes = 1:nSessions %Loop over all sessions:
    plot(params.xtime,sesmat(iSes,:,2,1),'r','LineWidth',0.25)
    plot(params.xtime,sesmat(iSes,:,2,2),'b','LineWidth',0.25)
end
%Figure make up:
title('Auditory Stim')
legend(params.areas)
xlim([-0.5e6 1e6])


%% Fig 3b - average population latency for each area:

%parameters for zoom in:
xmin = -0.02e6;
xmax = 0.08e6;
xresol = 20e3;
ymin = -1;
ymax = 3;
yresol = 1;

figure; hold all; set(gcf,'units','normalized','Position',[0.12 0.3 0.35 0.44],'color','w');

%Visual cortex: 
iArea = 2;
subplot(2,2,3); hold all;
%Plot the lines: 
meantoplot  = squeeze(nanmean(sesmat(:,:,:,iArea),1));
errortoplot = squeeze(nanstd(sesmat(:,:,:,iArea),[],1)) / sqrt(nSessions);
h = shadedErrorBar(params.xtime,meantoplot(:,2),errortoplot(:,2),{'r','LineWidth',0.5},0); delete(h.edge(:));
handles(1) = h.mainLine;
h = shadedErrorBar(params.xtime,meantoplot(:,1),errortoplot(:,1),{'b','LineWidth',0.5},0); delete(h.edge(:));
handles(2) = h.mainLine;
plot([xmin xmax xmax xmin xmin],[ymin ymin ymax ymax ymin],'k','LineWidth',0.5)

%Figure make up:
title('V1')
legend(handles,{'A' 'V'},'Location','NorthWest'); legend boxoff;
set(gca,'XTick',(-1e6:0.25e6:2e6),'XTickLabels',(-1e6:0.25e6:2e6)*1e-3)
xlim([-0.25e6 0.75e6])
ylim([-1 10])
xlabel('Time (ms)')
ylabel('Firing rate (z-score)')

%Close up figure of the rise of the population activity in V1:
subplot(2,2,1); hold all;
h = shadedErrorBar(params.xtime,meantoplot(:,2),errortoplot(:,2),{'r','LineWidth',0.5},0); delete(h.edge(:));
h = shadedErrorBar(params.xtime,meantoplot(:,1),errortoplot(:,1),{'b','LineWidth',0.5},0); delete(h.edge(:));
plot([xmin xmax],[1 1],':','Color',[0.4 0.4 0.4],'LineWidth',0.25)
xlim([xmin xmax])
ylim([ymin ymax])
set(gca,'XTick',xmin:xresol:xmax,'XTickLabels',(xmin:xresol:xmax)*1e-3,'YTick',ymin:yresol:ymax,'YTickLabels',(ymin:yresol:ymax),'XTickLabelRotation',45)

plot([poplat(iArea,1) poplat(iArea,1)]*1e3,[0 1],':','Color',[0 0 1],'LineWidth',0.25)
plot([poplat(iArea,2) poplat(iArea,2)]*1e3,[0 1],':','Color',[1 0 0],'LineWidth',0.25)
text((poplat(iArea,1)-5)*1e3,-0.1,[num2str(poplat(iArea,1)) ' ms'],'FontSize',12)
text((poplat(iArea,2)-5)*1e3,-0.3,[num2str(poplat(iArea,2)) ' ms'],'FontSize',12)

%Auditory cortex:
iArea = 1;
subplot(2,2,4); hold all;
meantoplot  = squeeze(nanmean(sesmat(:,:,:,iArea),1));
errortoplot = squeeze(nanstd(sesmat(:,:,:,iArea),[],1)) / sqrt(nSessions);
h = shadedErrorBar(params.xtime,meantoplot(:,2),errortoplot(:,2),{'r','LineWidth',0.5},0); delete(h.edge(:));
handles(1) = h.mainLine;
h = shadedErrorBar(params.xtime,meantoplot(:,1),errortoplot(:,1),{'b','LineWidth',0.5},0); delete(h.edge(:));
handles(2) = h.mainLine;
plot([xmin xmax xmax xmin xmin],[ymin ymin ymax ymax ymin],'k','LineWidth',0.5)

%Figure make up:
title('AC')
legend(handles,{'A' 'V'},'Location','NorthWest'); legend boxoff;
set(gca,'XTick',(-1e6:0.25e6:2e6),'XTickLabels',(-1e6:0.25e6:2e6)*1e-3)
xlim([-0.25e6 0.75e6])
ylim([-1 10])
xlabel('Time (ms)')
ylabel('Firing rate (z-score)')

subplot(2,2,2); hold all; handles = [];
h = shadedErrorBar(params.xtime,meantoplot(:,2),errortoplot(:,2),{'r','LineWidth',0.5},0); delete(h.edge(:));
h = shadedErrorBar(params.xtime,meantoplot(:,1),errortoplot(:,1),{'b','LineWidth',0.5},0); delete(h.edge(:));
plot([xmin xmax],[1 1],':','Color',[0.4 0.4 0.4],'LineWidth',0.25)
xlim([xmin xmax])
ylim([ymin ymax])
set(gca,'XTick',xmin:xresol:xmax,'XTickLabels',(xmin:xresol:xmax)*1e-3,'YTick',ymin:yresol:ymax,'YTickLabels',(ymin:yresol:ymax),'XTickLabelRotation',45)

plot([poplat(iArea,1) poplat(iArea,1)]*1e3,[0 1],':','Color',[0 0 1],'LineWidth',0.25)
plot([poplat(iArea,2) poplat(iArea,2)]*1e3,[0 1],':','Color',[1 0 0],'LineWidth',0.25)
text((poplat(iArea,1)-5)*1e3,-0.1,[num2str(poplat(iArea,1)) ' ms'],'FontSize',12)
text((poplat(iArea,2)-5)*1e3,-0.3,[num2str(poplat(iArea,2)) ' ms'],'FontSize',12)

filename = sprintf('PopulationLatency_A1V1_AV_AllCohorts_perarea.eps');
% export_fig(fullfile(params.savedir,filename),gcf);

%% Bootstrap test to compare auditory response latencies between AC and V1:
%Compute latency to cross z-scored threshold, averaged across sessions:
nBoot = 1000;
poplatboot = NaN(2,2,nBoot); %2x2 areas by modalities
for iArea = 1:2
    for iMod = 1:2
        for iBoot = 1:nBoot
            idx = randi(nSessions,nSessions,1);
            poplatboot(iArea,iMod,iBoot) = params.xtime(find(squeeze(nanmean(sesmat(idx,:,iMod,iArea),1))>params.zthreshold,1))*1e-3;
        end
    end
end

poplatboot(poplatboot<0) = NaN;
% Is the 97.5 slowest percentile of the bootstrap of AC auditory RT faster than the 2.5 percentile of V1?
hypo = prctile(squeeze(poplatboot(1,2,:)),97.5) < prctile(squeeze(poplatboot(2,2,:)),2.5) %percentile difference:
%This is not the case...

arealabels      = repmat({'AC', 'V1'}',nBoot*2,1);
modlabels       = repmat({'V', 'V','A','A'}',nBoot,1);
bootlabels      = repmat([1:1000],4,1);
bootlabels      = [bootlabels(:)];

tbl = table(poplatboot(:),arealabels,modlabels,bootlabels,...
    'VariableNames',{'Latency','Area','Modality','Bootiteration'}); %Create table for mixed model
writetable(tbl,'SourceData_Fig3b_AV_PopLatency.xlsx')

