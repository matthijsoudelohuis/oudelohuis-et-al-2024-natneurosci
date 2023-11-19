%% Oude Lohuis et al. 2023 Nat Neurosci
% Triple dissociation of auditory, visual, and motor processing in primary visual cortex
% MOL (C) 2023

% MOL_Conflict_RT
% Script analyzes the reaction times of the conflict trials

%% Paramters
params.savedir              = 'E:\Documents\PhD\Figures\Project CHDET\Results - auV1\8ConflictBehavior';

%% data
[Data] = MOL_GetData('E:','CHDET',{'BehaviorConflict'},{},[],{'sessionData' 'trialData'});
sessionData     = Data.sessionData;
trialData       = Data.trialData;

%% Remove last 20 trials:
trialData = MOL_RemoveLastnTrials(trialData,20);

%% Filter out sessions that do not have 4 levels of visual and auditory saliency:
idx = any([cellfun(@numel,sessionData.vecFreqChange) cellfun(@numel,sessionData.vecOctChange)]==5,2) & ...
any(cellfun(@numel,sessionData.vecOriChange)==4,2);
[sessionData, trialData]        = MOL_getTempPerSes(sessionData.session_ID(idx),sessionData,trialData);

%% Save dataset:
save('Dataset6_2.mat','params','sessionData','trialData')

%% Or load dataset
load Dataset6_2.mat


%% Make figure of response latencies to different trial types with conflict trials:
splits          = {};
params.trialcolors = {};

for iSplit = 2:5
    splits{end+1}               = strcmp(trialData.trialType,'Y') & trialData.audioFreqChangeNorm==iSplit & trialData.vecResponse==1 & ~(trialData.hasphotostim==1);
    params.trialcolors{end+1}   = [(iSplit-1)/4 0 0];
end

for iSplit = 2:5
    splits{end+1}               = strcmp(trialData.trialType,'X') & trialData.visualOriChangeNorm==iSplit & trialData.vecResponse==2 & ~(trialData.hasphotostim==1);
    params.trialcolors{end+1}   = [0 0 (iSplit-1)/4];
end

for iSplit = 2:5
    splits{end+1}               = strcmp(trialData.trialType,'C') & trialData.audioFreqChangeNorm==iSplit & trialData.visualOriChangeNorm==iSplit & trialData.vecResponse==1 & ~(trialData.hasphotostim==1);
    params.trialcolors{end+1}   = [(iSplit-1)/4 0 (iSplit-1)/8];
end

for iSplit = 2:5
    splits{end+1}               = strcmp(trialData.trialType,'C') & trialData.audioFreqChangeNorm==iSplit & trialData.visualOriChangeNorm==iSplit & trialData.vecResponse==2 & ~(trialData.hasphotostim==1);
    params.trialcolors{end+1}   = [(iSplit-1)/8 0 (iSplit-1)/4];
end

% params.triallabels = {'Asub' 'Athr' 'Asup' 'Amax' 'Vsub' 'Vthr' 'Vsup' 'Vmax' 'CsubA' 'CsubV' 'CthrA' 'CthrV' 'CsupA' 'CsupV' 'CmaxA' 'CmaxV'};
params.triallabels = repmat({'Sub' 'Thr' 'Sup' 'Max'},1,4);

nSplits         = length(splits);
meantoplot      = NaN(nSplits,1);
errortoplot      = NaN(nSplits,1);

for iSplit = 1:nSplits
    meantoplot(iSplit)    = nanmean(trialData.responseLatency(splits{iSplit}));
    errortoplot(iSplit)   = nanstd(trialData.responseLatency(splits{iSplit}));
    errortoplot(iSplit)   = errortoplot(iSplit) / sqrt(sum(splits{iSplit}));
end

meantoplot = meantoplot*1e-3; %convert to ms
errortoplot = errortoplot*1e-3; %convert to ms

figure; hold all; set(gcf,'units','normalized','Position',[0.1 0.25 0.23 0.27],'color','w');
xlim([0.5 nSplits+0.5])
for iSplit = 1:length(splits)
    z = errorbar(iSplit,meantoplot(iSplit),errortoplot(iSplit),'.','Color','k','MarkerSize',25,'MarkerFaceColor',params.trialcolors{iSplit},'MarkerEdgeColor',params.trialcolors{iSplit},'LineWidth',1.5,'CapSize',0);
end
set(gca,'XTick',1:nSplits,'XTickLabel',params.triallabels,'XTickLabelRotation',45,'FontSize',10)
set(gca,'YTick',[300 425 550],'FontSize',10)

handles = [];
for iSplit = 1:4
    handles(iSplit) = errorbar(25,meantoplot(iSplit),0,'.','Color','k','MarkerSize',35,'MarkerFaceColor',[iSplit/5 iSplit/5 iSplit/5],'MarkerEdgeColor',[iSplit/5 iSplit/5 iSplit/5],'LineWidth',1,'CapSize',0);
end
legend(handles,{'Sub' 'Thr' 'Sup' 'Max'},'Location','NorthWest','FontSize',12); legend boxoff;

ylabel('Reaction time (ms)','FontSize',10)
ylim([300 550])
filename = sprintf('ConflictRT.eps');
export_fig(fullfile(params.savedir,filename),gcf)

%% Difference in RT between Au and Vis of saliency matched conditions: 

rtdiff = nanmean(meantoplot(5:8) - meantoplot(1:4));
fprintf('The difference in reaction time between saliency matched auditory and visual conditions is on average: %3.1f ms\n',rtdiff)

