startover

%% Parameter settings
params.Experiments          = {'ChangeDetectionConflictDecor' 'VisOnlyTwolevels' 'ChangeDetectionConflict' }; %Which versions of the task to load data from
%% Get input arguments:
[Data] = MOL_GetData('E:','CHDET',params.Experiments,{},{},{'sessionData' 'trialData_newtrials' 'spikeData' 'videoData'});
sessionData     = Data.sessionData;
trialData       = Data.trialData;
spikeData       = Data.spikeData;
videoData       = Data.videoData;

%% Filter out sessions with muscimol:
sesids              = sessionData.session_ID(cellfun(@isempty,sessionData.MuscimolArea));
fprintf('Removed %d/%d sessions  with muscimol manipulations\n',length(sessionData.session_ID)-length(sesids),length(sessionData.session_ID));
[sessionData, trialData, spikeData,videoData]        = MOL_getTempPerSes(sesids,sessionData,trialData,spikeData,videoData);

%% Filter out sessions with no systematic post orientation or frequency:
sesids = unique(trialData.session_ID(~isnan(trialData.visualOriPostChangeNorm) & ~isnan(trialData.audioFreqPostChangeNorm)));
[sessionData, trialData, spikeData,videoData]        = MOL_getTempPerSes(sesids,sessionData,trialData,spikeData,videoData);

%% Only sessions with video:
% Filter sessions based on containing videoData:
[sessionData,trialData,spikeData,videoData] = MOL_getTempPerSes(unique(videoData.session_ID),sessionData,trialData,spikeData,videoData);

% Filter sessions based on containing motSVD data:
[sessionData,trialData,spikeData,videoData] = MOL_getTempPerSes(videoData.session_ID(~cellfun(@isempty,videoData.motSVD)),sessionData,trialData,spikeData,videoData);

%% Filter out sessions that are active:
sesids = unique(sessionData.session_ID(strcmp(sessionData.State,'Behaving')));
[sessionData, trialData, spikeData,videoData]        = MOL_getTempPerSes(sesids,sessionData,trialData,spikeData,videoData);

%% Report dataset:
fprintf('Dataset: %d sessions, %d trials, %d neurons, %d videos\n',length(sessionData.session_ID),length(trialData.session_ID),length(spikeData.session_ID),length(videoData.session_ID));
nSessions = length(sessionData.session_ID);

%%
sessionData.Rec_datetime
