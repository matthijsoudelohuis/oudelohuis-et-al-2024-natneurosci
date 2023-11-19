%% Script that shows spiking activity and video motion following auditory stimuli
% in a audiovisual change detection task
% MOL (C) 2023
% as reported in Oude Lohuis et al. 2023 Nat Neurosci
% "Triple dissociation of auditory, visual and motor processing in primary visual cortex"

startover

%% Parameter settings:
params                      = params_histresponse_auV1;
params.Experiments          = {'ChangeDetectionConflictDecor' 'VisOnlyTwolevels' 'ChangeDetectionConflict' }; %Which versions of the task to load data from
params.ExperimentLabels     = {'NE' 'UST' 'MST'}; %Which versions of the task to load data from
params.AlignOn              = 'stimChange';         %On which timestamp to align trials to
params.nExperiments         = length(params.Experiments);
params                      = MOL_getColors_CHDET(params);

params.videofield           = 'motSVD';

params.areas                = {'V1'};
params.nAreas               = length(params.areas);

params.savedir              = 'E:\Documents\PhD\Figures\Project CHDET\Results - auV1\6FeatureCoding\ExampleNeurons\';

%% Get data:
% [Data] = MOL_GetData('E:','CHDET',params.Experiments,{'1008' '1009' '2003' '2009' '2010' '2011' '2012' '2030' '2031'},[],{'sessionData' 'trialData_newtrials' 'spikeData' 'videoData'});
[Data] = MOL_GetData('E:','CHDET',params.Experiments,{'2009' '2010' '2011' '2012' '2030' '2031'},[],{'sessionData' 'trialData_newtrials' 'spikeData' 'videoData'});
sessionData         = Data.sessionData;
trialData           = Data.trialData;
spikeData           = Data.spikeData;
videoData           = Data.videoData;

%% Remove last 20 trials:
trialData           = MOL_RemoveLastnTrials(trialData,20);

%% Filter out neurons based on quality:
spikeData           = MOL_filterNeurons(sessionData,trialData,spikeData);

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
[sessionData,trialData,spikeData,videoData] = MOL_getTempPerSes(videoData.session_ID(~cellfun(@isempty,videoData.motSVD)),sessionData,trialData,spikeData,videoData);

%% Report dataset:
fprintf('Dataset: %d sessions, %d trials, %d neurons, %d videos\n',length(sessionData.session_ID),length(trialData.session_ID),length(spikeData.session_ID),length(videoData.session_ID));
nSessions = length(sessionData.session_ID);

%% Save dataset:
save('Dataset2_2.mat','params','sessionData','trialData')

%% Or load dataset:
load Dataset2_2.mat



%% Parameters for video hist:
params.t_pre        = -1e6;
params.t_post       = 4e6;
params.nSVDs        = 25;

[video_hist,video_hist_tot,video_hist_z,params] = MOL_PSTH_video(sessionData,trialData,videoData,params);

%% Parameters for raster plots:
params.exportfig    = 0;
params.t_pre        = -0.4e6;
params.t_post       = 0.8e6;
params.zscore       = 0;
params.minNtrials   = 5;

params.smoothSVD    = 3;

%% Show raster plot of example neurons showing correlation to movements irrespective of saliency:

cell_IDs            = {};
cell_IDs{end+1}     = '20092018082331113'; %Paper example #1

%Other examples:
% cell_IDs{end+1}     = '10092019031211313';
% cell_IDs{end+1}     = '10092019030831126';

params.trialcategories = 'Saliency';

MOL_plotRasterVid_auV1(sessionData,trialData,spikeData,video_hist_z,cell_IDs,params)


%% Show raster plot of example neurons showing auditory responses irrespective of movement:

cell_IDs            = {};
% cell_IDs{end+1}     = '20122018081611121'; %former paper example: unrelated to movement
cell_IDs{end+1}     =     '20102018081611378'; %new paper example

% cell_IDs{end+1}     = '10122019041031329'; %Auditory no motor example cell glm
 
% %Other examples:
% cell_IDs{end+1}     = '20102018081611325';
% cell_IDs{end+1}     = '20122018081431167'; %fast transient
% cell_IDs{end+1}     = '20112018081011160'; %fast transient
% cell_IDs{end+1}     = '20092018082211178'; %sharp onset, not related to total movement
% cell_IDs{end+1}     = '20092018082211082'; %tuned sharp onset, not related to total movement
% cell_IDs{end+1}     = '10092019031211324';
% cell_IDs{end+1}     = '10092019031211313'; %sharp au onset not relate to movement
% cell_IDs{end+1}     = '10092019031311154'; %fast transient
% cell_IDs{end+1}     = '10082019031211237'; %sharp au onset not relate to movement
% cell_IDs{end+1}     = '10092019030831126';
% cell_IDs{end+1}     = '20122018081431167'; %fast transient
% cell_IDs{end+1}     = '20092018082411239'; % %Very sharp onset transet
% cell_IDs{end+1}     = '20122018081431160'; %  %tuned?
% cell_IDs{end+1}     = '20102018081311171'; %

params.trialcategories = 'Saliency';
params.exportfig = 1;
MOL_plotRasterVid_auV1(sessionData,trialData,spikeData,video_hist_z,cell_IDs,params)


%% Show raster plot of example neurons showing orientation and frequency selective responses:
cell_IDs            = {};

%Orientation selective example neuron:
cell_IDs{end+1}     = '20302020011621372'; %Paper example: primary example with persistent orientation selectivity

%Other examples:
% cell_IDs{end+1}     = '20302020011621304';
% cell_IDs{end+1}     = '20312020012311114';
% cell_IDs{end+1}     = '20312020012311127';
% cell_IDs{end+1}     = '20312020012311256';

% Frequency selective responses:
cell_IDs{end+1}     = '20122018081431146'; %paper example: strong tuning
cell_IDs{end+1}     = '20112018081011160'; %paper example: tuning related to movements

%Other examples:
% cell_IDs{end+1}     = '20092018082211082';
% cell_IDs{end+1}     = '20092018082331113';
% cell_IDs{end+1}     = '20092018082411084';
% cell_IDs{end+1}     = '20122018081311186';
% cell_IDs{end+1}     = '20122018081431040';
% cell_IDs{end+1}     = '20122018081431094';
% cell_IDs{end+1}     = '20122018081431095';
% cell_IDs{end+1}     = '20122018081431167'; 
% cell_IDs{end+1}     = '20122018081431176';

params.trialcategories = 'PostFeature2';

MOL_plotRasterVid_auV1(sessionData,trialData,spikeData,video_hist_z,cell_IDs,params)




