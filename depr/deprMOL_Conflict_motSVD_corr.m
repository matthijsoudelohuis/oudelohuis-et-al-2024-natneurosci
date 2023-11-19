startover

%% Parameter settings
params                      = params_histresponse_auV1;% All time is in microseconds

params.Experiments          = {'ChangeDetectionConflict' }; %Which versions of the task to load data from
params.ExperimentLabels     = {'MST'}; %Labels for the different experiments
params.nExperiments         = length(params.Experiments);

params.AlignOn              = 'stimChange';      %On which timestamp to align as t=0

params.videofield           = 'motSVD';

params.area                 = 'V1';

params.fs                   = 25; %Hz

params.t_pre                = -1e6;
params.t_post               = 4e6;

params                      = MOL_getColors_CHDET(params);

params.savedir              = 'E:\Documents\PhD\Figures\Project CHDET\Results - auV1\8ConflictBehavior\ConflictMovements\';

%Set colors and labels for conditions:
params.colors_splits = [];
params.colors_splits(1,1:2,:) = [160 160 160; 143 43 44];% 157 0 214]; 
params.colors_splits(2,1:2,:) = [37 45 138; 157 0 214];% 157 0 214]; 
params.colors_splits = params.colors_splits / 256;
params.labels_splits = {'av' 'Av'; 'aV' 'AV'};
params.lines_splits = {'-'  '--'; '-.' ':'};

%% Get input arguments:
[Data] = MOL_GetData('E:','CHDET',params.Experiments,{},{},{'sessionData' 'trialData_newtrials' 'videoData'});
% [Data] = MOL_GetData('E:','CHDET',{'ChangeDetectionConflictDecor'},{'1008' '1009'},{},{'sessionData' 'trialData' 'spikeData' 'videoData'});
sessionData     = Data.sessionData;
trialData       = Data.trialData;
% spikeData       = Data.spikeData;
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

%% Filter out sessions with muscimol:
sesids              = sessionData.session_ID(cellfun(@isempty,sessionData.MuscimolArea));
fprintf('Removed %d/%d sessions  with muscimol manipulations\n',length(sessionData.session_ID)-length(sesids),length(sessionData.session_ID));
[sessionData, trialData, videoData]        = MOL_getTempPerSes(sesids,sessionData,trialData,videoData);

%% Filter out sessions with no systematic post orientation or frequency:
sesids = unique(trialData.session_ID(~isnan(trialData.visualOriPostChangeNorm) & ~isnan(trialData.audioFreqPostChangeNorm)));
[sessionData, trialData, videoData]        = MOL_getTempPerSes(sesids,sessionData,trialData,videoData);

%% Only sessions with video:
% Filter sessions based on containing videoData:
[sessionData,trialData,videoData] = MOL_getTempPerSes(unique(videoData.session_ID),sessionData,trialData,videoData);

% Filter sessions based on containing motSVD data:
[sessionData,trialData,videoData] = MOL_getTempPerSes(videoData.session_ID(~cellfun(@isempty,videoData.motSVD)),sessionData,trialData,videoData);

%% Filter sessions that have one or four postFreq and postOri during conflict trials:
nSessions           = length(sessionData.session_ID);
singleStimConfl     = false(nSessions,1);
for iSes = 1:nSessions
    postfreqs               = trialData.audioFreqPostChangeNorm(strcmp(trialData.trialType,'C') & strcmp(trialData.session_ID,sessionData.session_ID(iSes)));
    postoris                = trialData.visualOriPostChangeNorm(strcmp(trialData.trialType,'C') & strcmp(trialData.session_ID,sessionData.session_ID(iSes)));
    singleStimConfl(iSes)   = ismember(numel(unique(postfreqs)),[1 4]) && ismember(numel(unique(postoris)),[1 4]);
%     singleStimConfl(iSes)   = ismember(numel(unique(postfreqs)),1) && ismember(numel(unique(postoris)),1);
end
fprintf('Removed %d/%d sessions with multiple conflict stimuli\n',length(sessionData.session_ID)-sum(singleStimConfl),length(sessionData.session_ID));
[sessionData,trialData,videoData] = MOL_getTempPerSes(sessionData.session_ID(singleStimConfl),sessionData,trialData,videoData);

%% Report dataset:
fprintf('Dataset: %d sessions, %d trials, %d videos \n',length(sessionData.session_ID),length(trialData.session_ID),length(videoData.session_ID));
nSessions = length(sessionData.session_ID);

%%




%adapt to use the mol psth video function !!!!



%% Shift all motion SVD values given that motion is computed as diff between frames
for iSes = 1:nSessions
    videoData.motSVD{iSes} = circshift(videoData.motSVD{iSes},2,1);
end

%%
params.t_pre        = -1e6;
params.t_post       = 4e6;
params.nSVDs        = 25;

[video_hist,video_hist_tot,video_hist_z,params] = MOL_PSTH_video(sessionData,trialData,videoData,params);


% %% Shift all motion SVD values given that motion is computed as diff between frames
% fprintf('Shifting frames for session    \n');
% for iSes = 1:nSessions
%     fprintf(repmat('\b', 1, numel([num2str(iSes-1) num2str(nSessions)])+1));
%     fprintf('%d/%d',iSes,nSessions);
% %     videoData.motSVD{iSes} = [nan(1,500); videoData.motSVD{iSes}(1:end-1,:)];
%     videoData.motSVD{iSes} = circshift(videoData.motSVD{iSes},2,1);
% end
% fprintf('\n')
% 
% %% Set variables for video processing:
% params.videofield       = 'motSVD';
% 
% edges                   = params.t_pre:1e6/params.fs:params.t_post;
% params.xtime_video      = edges+1e6/params.fs/2;
% params.nTimebins_video  = numel(edges);
% nTrials             	= length(trialData.stimChange);


% %% Make peri-stimulus histogram of motion
% params.videofield = 'motSVD';
% 
% edges                       = params.t_pre:1e6/params.fs:params.t_post;
% % params.xtime    = edges(1:end-1)+1e6/params.fs/2;
% params.xtime_video          = edges+1e6/params.fs/2;
% params.nTimebins_video      = numel(edges);
% params.nSVDs                = 25;
% 
% hist_mat                    = NaN(params.nTimebins_video,params.nSVDs,nTrials);
% 
% fprintf('Computing svd response for trial        \n');
% for iTrial = 1:nTrials
%     fprintf(repmat('\b', 1, numel([num2str(iTrial-1) num2str(nTrials)])+1));
%     fprintf('%d/%d',iTrial,nTrials);
%     ses_idx                             = strcmp(sessionData.session_ID,trialData.session_ID(iTrial));
%     
%     idx                                 = videoData.ts{ses_idx}>trialData.(params.AlignOn)(iTrial)+params.t_pre & videoData.ts{ses_idx}<trialData.(params.AlignOn)(iTrial)+params.t_post;
%     hist_mat(:,:,iTrial)                = interp1(videoData.ts{ses_idx}-trialData.(params.AlignOn)(iTrial),videoData.(params.videofield){ses_idx}(:,1:params.nSVDs),params.xtime_video,'linear');
% end
% fprintf('\n')



%%
%Parameters correlation:
params.twin_resp_start  = 0e6;
params.twin_resp_stop   = 2e6;
idx_time                = params.xtime_video>params.twin_resp_start & params.xtime_video<=params.twin_resp_stop; %index across timebins across which correlation of termporal firing rate will be calculated

nVsals                  = 2; %number of visual saliencies (thr and max)
nAsals                  = 2; %number of auditory saliencies (thr and max)

Rvs                     = NaN(nSessions,nVsals,nAsals); %init variables (Rvs is correlation of movements during conflict trial to stimulus-matched visual unimodal trial)
Ras                     = NaN(nSessions,nVsals,nAsals);

% Rvsp                    = NaN(nSessions,nVsals,nAsals);
% Rasp                    = NaN(nSessions,nVsals,nAsals);
 params.minNtrialCond = 5;
 
fprintf('Computing correlation of conflict trials to unimodal trial types           \n');

for iSes = 1:nSessions %Loop over all sessions:
    fprintf(repmat('\b', 1, numel([num2str(iSes-1) num2str(nSessions)])+2));
    fprintf('%d/%d\n',iSes,nSessions);
    
    %Get the relevant data for each session individually:
    temptrialData               = MOL_getTempPerSes(sessionData.session_ID(iSes),trialData);
    idx_ses                     = strcmp(trialData.session_ID,sessionData.session_ID(iSes));
    
    %subselect motion psthistogram:
    hist_mat_temp               = hist_mat(:,:,idx_ses);
    
    %See if session has one stimulus condition for conflict trials or multiple:
    if numel(unique(temptrialData.visualOriPostChangeNorm(strcmp(temptrialData.trialType,'C'))))==1
        uniqueflag = true;
    elseif numel(unique(temptrialData.visualOriPostChangeNorm(strcmp(temptrialData.trialType,'C'))))==4
        uniqueflag = false;
        temptrialData.visualOriPostChangeNorm = round(temptrialData.visualOriPostChangeNorm/2);
        temptrialData.audioFreqPostChangeNorm = round(temptrialData.audioFreqPostChangeNorm/2);
    else
        error('unknown number of postchange features during conflict trials')
    end
    
    for iVsal = 1:nVsals
        for iAsal = 1:nAsals

            if uniqueflag
                idx_C   = strcmp(temptrialData.trialType,'C') & temptrialData.visualOriChangeNorm==iVsal+1 & temptrialData.audioFreqChangeNorm==iAsal+1;
                idx_V   = strcmp(temptrialData.trialType,'X') & temptrialData.visualOriChangeNorm==iVsal+1 & temptrialData.visualOriPostChangeNorm==unique(temptrialData.visualOriPostChangeNorm(idx_C));
                idx_A   = strcmp(temptrialData.trialType,'Y') & temptrialData.audioFreqChangeNorm==iAsal+1 & temptrialData.audioFreqPostChangeNorm==unique(temptrialData.audioFreqPostChangeNorm(idx_C));
                
                tempV = nanmean(hist_mat_temp(idx_time,:,idx_V),3);
                tempA = nanmean(hist_mat_temp(idx_time,:,idx_A),3);
                tempC = nanmean(hist_mat_temp(idx_time,:,idx_C),3);
                
                Rvs(iSes,iVsal,iAsal)                = corr(tempV(:),tempC(:),'rows','complete');
                Ras(iSes,iVsal,iAsal)                = corr(tempA(:),tempC(:),'rows','complete');
                
            else
                postclass = 1;
                idx_C1   = strcmp(temptrialData.trialType,'C') & temptrialData.visualOriChangeNorm==iVsal+1 & temptrialData.audioFreqChangeNorm==iAsal+1 & ...
                    temptrialData.visualOriPostChangeNorm==postclass & temptrialData.audioFreqPostChangeNorm==postclass;
                idx_V1   = strcmp(temptrialData.trialType,'X') & temptrialData.visualOriChangeNorm==iVsal+1 & temptrialData.visualOriPostChangeNorm==postclass;
                idx_A1   = strcmp(temptrialData.trialType,'Y') & temptrialData.audioFreqChangeNorm==iAsal+1 & temptrialData.audioFreqPostChangeNorm==postclass;
                postclass = 2;
                idx_C2   = strcmp(temptrialData.trialType,'C') & temptrialData.visualOriChangeNorm==iVsal+1 & temptrialData.audioFreqChangeNorm==iAsal+1 & ...
                    temptrialData.visualOriPostChangeNorm==postclass & temptrialData.audioFreqPostChangeNorm==postclass;
                idx_V2   = strcmp(temptrialData.trialType,'X') & temptrialData.visualOriChangeNorm==iVsal+1 & temptrialData.visualOriPostChangeNorm==postclass;
                idx_A2   = strcmp(temptrialData.trialType,'Y') & temptrialData.audioFreqChangeNorm==iAsal+1 & temptrialData.audioFreqPostChangeNorm==postclass;
                
                tempV1 = nanmean(hist_mat_temp(idx_time,:,idx_V1),3);
                tempA1 = nanmean(hist_mat_temp(idx_time,:,idx_A1),3);
                tempC1 = nanmean(hist_mat_temp(idx_time,:,idx_C1),3);
                
                tempV2 = nanmean(hist_mat_temp(idx_time,:,idx_V2),3);
                tempA2 = nanmean(hist_mat_temp(idx_time,:,idx_A2),3);
                tempC2 = nanmean(hist_mat_temp(idx_time,:,idx_C2),3);
                
                Rvs(iSes,iVsal,iAsal)                = corr([tempV1(:); tempV2(:)],[tempC1(:); tempC2(:)],'rows','complete');
                Ras(iSes,iVsal,iAsal)                = corr([tempA1(:); tempA2(:)],[tempC1(:); tempC2(:)],'rows','complete');
                
            end
        end
    end
end


%% Figure of histogram of correlations to visual and auditory stim components / decision components:
figure; hold all; set(gcf,'units','normalized','Position',[0.15 0.5 0.12 0.24],'color','w');

h = violinplot([Rvs(:) Ras(:)],[],'ViolinAlpha',1,'EdgeColor',[1 1 1],'BandWidth',0.15,'BoxColor',[0 0 0]);
h(1).ViolinColor = squeeze(params.colors_splits(2,1,:));
h(2).ViolinColor = squeeze(params.colors_splits(1,2,:));

ylim([-0.5 1])
p = signrank(Rvs(:),Ras(:));
sigstar([1 2],p)
set(gca,'XTick',[],'YTick',[-0.5 0 0.5 1])

fprintf('Correlation motion SVDs during conflict trials with unimodal trials types:\n')
fprintf('Significantly more correlated to auditory than visual unimodal trials \n (%1.2f vs %1.2f, p=%1.4f, Wilcoxon Signed Rank test) \n',nanmedian(Rvs(:)),nanmedian(Ras(:)),p)

filename = sprintf('motSVD_corrConflict.eps');
export_fig(fullfile(params.savedir,filename),gcf);
