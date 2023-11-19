%% Oude Lohuis et al. 2023 Nat Neurosci
% Triple dissociation of auditory, visual, and motor processing in primary visual cortex
% MOL (C) 2023

%% Get input arguments:
[Data] = MOL_GetData('E:','CHDET',{'BehaviorConflict'},{},{},{'sessionData' 'trialData' 'videoData'});
sessionData     = Data.sessionData;
trialData       = Data.trialData;
videoData       = Data.videoData;

%% Remove last 20 trials:
trialData = MOL_RemoveLastnTrials(trialData,20);

%% Filter out sessions that do not have 4 levels of visual and auditory saliency:
idx = any([cellfun(@numel,sessionData.vecFreqChange) cellfun(@numel,sessionData.vecOctChange)]==5,2) & ...
any(cellfun(@numel,sessionData.vecOriChange)==4,2);
[sessionData, trialData, videoData]        = MOL_getTempPerSes(sessionData.session_ID(idx),sessionData,trialData,videoData);

%% Filter out trials with photostimulation in V1:
sesids              = sessionData.session_ID(sessionData.UseOpto & (strcmp(sessionData.PhotostimArea,'V1') | strcmp(sessionData.PhotostimArea,'PPC')));
idx                 = trialData.hasphotostim==1 & ismember(trialData.session_ID,sesids);
trialFields         = fieldnames(trialData);
fprintf('Removed %d/%d trials  with photostim in V1\n',sum(idx),length(trialData.session_ID));
for iF = 1:length(trialFields)
    trialData.(trialFields{iF}) = trialData.(trialFields{iF})(~idx);
end

%% Filter sessions based on containing videoData:
[sessionData,trialData,videoData] = MOL_getTempPerSes(unique(videoData.session_ID),sessionData,trialData,videoData);

%% Save dataset:
save('Dataset6_3.mat','sessionData','trialData','videoData')

%% Or load dataset
load Dataset6_3.mat

%% Report dataset:
fprintf('Dataset: %d sessions, %d trials, %d videos\n',length(sessionData.session_ID),length(trialData.session_ID),length(videoData.session_ID));

%% Parameter settings for PSTH
params                      = params_histresponse(); % All time is in microseconds
params.histmethod           = 'individual'; %Whether to bin (and smooth) individual trials 'individual' or or all trials together ('total')
params.zscore               = 1;

params.SortBy               = 'maxResponse';
params.AlignOn              = 'stimChange';      %On which timestamp to align as t=0

params.videofield           = 'zarea';

params.fs                   = 25; %Hz

params.t_pre                = -1e6;
params.t_post               = 4e6;

params                      = MOL_getColors_CHDET(params);

% params.triallines           = {'-.' '-' '-.' '-' '-.' '-' '-.' '-'};
params.triallabels          = {'Sub' 'Thr' 'Sup' 'Max' 'Sub' 'Thr' 'Sup' 'Max'};


%% Compute peristimulus histogram of pupil size:
edges           = params.t_pre:1e6/params.fs:params.t_post;
% params.xtime    = edges(1:end-1)+1e6/params.fs/2;
params.xtime    = edges+1e6/params.fs/2;
nTimebins       = numel(edges);
nTrials         = length(trialData.stimChange);

hist_mat        = NaN(nTrials,nTimebins);
hist_mat_isgood = NaN(nTrials,nTimebins);

fprintf('Computing Z-scored response for trial        \n');
for iTrial = 1:nTrials
    fprintf(repmat('\b', 1, numel([num2str(iTrial-1) num2str(nTrials)])+1));
    fprintf('%d/%d',iTrial,nTrials);
    ses_idx                             = strcmp(sessionData.session_ID,trialData.session_ID(iTrial));
    idx                                 = videoData.ts{ses_idx}>trialData.(params.AlignOn)(iTrial)+params.t_pre & videoData.ts{ses_idx}<trialData.(params.AlignOn)(iTrial)+params.t_post;
    hist_mat(iTrial,1:sum(idx))     	= videoData.zarea{ses_idx}(idx);
    %     hist_mat(iTrial,1:sum(idx))     	= videoData.area{ses_idx}(idx);
    hist_mat_isgood(iTrial,1:sum(idx))  = videoData.isgood{ses_idx}(idx);
end
fprintf('\n')

%% Fig. S10a: mean pupil response to the different saliencies:
params.trialcolors = {};

splits          = {};
for iSplit = 2:5
    splits{end+1}               = strcmp(trialData.trialType,'Y') & trialData.audioFreqChangeNorm==iSplit;
    params.trialcolors{end+1}   = [iSplit/5 0 0];
end

for iSplit = 2:5
    splits{end+1}               = strcmp(trialData.trialType,'X') & trialData.visualOriChangeNorm==iSplit;
    params.trialcolors{end+1}   = [0 0 iSplit/5];
end

meantoplot = NaN(length(splits),nTimebins);
errortoplot = NaN(length(splits),nTimebins);

for iSplit = 1:length(splits)
    meantoplot(iSplit,:)    = nanmean(hist_mat(splits{iSplit},:),1);
    meantoplot(iSplit,:)    = meantoplot(iSplit,:) - mean(meantoplot(iSplit,params.xtime<0));
    
    errortoplot(iSplit,:)   = nanstd(hist_mat(splits{iSplit},:),1);
    errortoplot(iSplit,:)   = errortoplot(iSplit,:) / sqrt(sum(splits{iSplit}));
end

handles = [];
figure; hold all; set(gcf,'units','normalized','Position',[0.05 0.4 0.36 0.21],'color','w');
subplot(1,2,1)
for iSplit = 1:4%1:length(splits)
    h = shadedErrorBar(params.xtime*1e-6,meantoplot(iSplit,:),errortoplot(iSplit,:),{'-','markerfacecolor',params.trialcolors{iSplit},'LineWidth',3},0);
    h.mainLine.Color = params.trialcolors{iSplit};    h.patch.FaceColor = params.trialcolors{iSplit};
    delete(h.edge(1)); delete(h.edge(2));
    handles(end+1) = h.mainLine; hold all;
end
xlabel('Time (in sec)')
ylabel('Z-scored pupil area')
xlim([edges(1)*1e-6 edges(end-1)*1e-6])
ylim([-0.2 1])
legend(fliplr(handles),params.triallabels(4:-1:1),'Location','NorthWest','FontSize',10); legend boxoff;

handles = [];
subplot(1,2,2)
for iSplit = 5:8%1:length(splits)
    h = shadedErrorBar(params.xtime*1e-6,meantoplot(iSplit,:),errortoplot(iSplit,:),{'-','markerfacecolor',params.trialcolors{iSplit},'LineWidth',3},0);
    h.mainLine.Color = params.trialcolors{iSplit};    h.patch.FaceColor = params.trialcolors{iSplit};
    delete(h.edge(1)); delete(h.edge(2));
    handles(end+1) = h.mainLine; hold all;
end

xlabel('Time (in sec)')
ylabel('Z-scored pupil area')
xlim([edges(1)*1e-6 edges(end-1)*1e-6])
ylim([-0.2 1])
legend(fliplr(handles),params.triallabels(8:-1:5),'Location','NorthWest','FontSize',10); legend boxoff;
MOL_prepfigAI

%%
params.t_respwin_start  = 0;
params.t_respwin_stop   = 3e6;

hist_mat_baselinecorr       = hist_mat - repmat(mean(hist_mat(:,params.xtime<0),2),1,nTimebins);
maxresp                     = max(hist_mat_baselinecorr(:,params.xtime>params.t_respwin_start & params.xtime<params.t_respwin_stop),[],2);

meantoplot                  = NaN(length(splits),1);
errortoplot                 = NaN(length(splits),1);

for iSplit = 1:length(splits)
    meantoplot(iSplit)    = nanmean(maxresp(splits{iSplit},:),1);
    errortoplot(iSplit)   = nanstd(maxresp(splits{iSplit},:),1);
    errortoplot(iSplit)   = errortoplot(iSplit,:) / sqrt(sum(splits{iSplit}));
end

figure; hold all; set(gcf,'units','normalized','Position',[0.45 0.35 0.17 0.24],'color','w');
h = [];
for iSplit = 1:length(splits)
    h(iSplit) = bar(iSplit,meantoplot(iSplit),0.8);
    set(h(iSplit),'facecolor',params.trialcolors{iSplit});
    z = errorbar(iSplit,meantoplot(iSplit),errortoplot(iSplit),'k.','LineWidth',2);
    errorbar_tick(z,0.001,'units')
end

set(gca,'XTick',1:8,'XTickLabel',params.triallabels,'XTickLabelRotation',45)
set(gca,'YTick',[0 0.5  1 1.5])
ylabel('Pupil dilation (z-score)')
ylim([0 1.5])

%% Statistics:

G_mou           = cell(nTrials,1);
uMice           = unique(sessionData.mousename);
for iMouse = 1:length(uMice)
    G_mou(ismember(trialData.session_ID,sessionData.session_ID(strcmp(sessionData.mousename,uMice{iMouse})))) = uMice(iMouse);
end

idx             = ismember(trialData.trialType,{'X' 'Y'});

X_sal           = [trialData.visualOriChangeNorm trialData.audioFreqChangeNorm];
X_sal           = X_sal-1;
X_sal           = sum(X_sal,2);

X_mod           = [trialData.hasvisualchange trialData.hasaudiochange*2];
X_mod           = sum(X_mod,2);

X_hit           = [trialData.correctResponse];

tbl             = table(maxresp(idx),X_mod(idx),X_sal(idx),X_hit(idx),G_mou(idx),'VariableNames',{'Pupil','Modality','Saliency','Hit','Mouse'}); %Create table for mixed model

lme             = fitlme(tbl,'Pupil~Modality+Saliency+Hit+(1|Mouse)'); %construct linear mixed effects model with fixed effect of temporal window (NO random intercept for mice as cohorts are tested)
stats           = dataset2table(anova(lme,'DFMethod','Satterthwaite')); %Perform ANOVA on model and output as matrix
fprintf('Effect modality on pupil dilation: \n')
fprintf('F(%d,%2.0f)=%1.2f, p=%1.4f; \n',stats{2,3},stats{2,4},stats{2,2},stats{2,5})
fprintf('Effect saliency on pupil dilation: \n')
fprintf('F(%d,%2.0f)=%1.2f, p=%1.2e; \n',stats{3,3},stats{3,4},stats{3,2},stats{3,5})
fprintf('Effect hit/miss on pupil dilation: \n')
fprintf('F(%d,%2.0f)=%1.2f, p=%1.2e; \n',stats{4,3},stats{4,4},stats{4,2},stats{4,5})

writetable(tbl,'SourceData_FigS10b_PupilDilation.xlsx')


