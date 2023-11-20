%% Script that shows spiking activity and video motion following auditory stimuli
% in a audiovisual change detection task
% MOL (C) 2023
% Oude Lohuis et al. 2024 Nat Neurosci
% "Triple dissociation of auditory, visual and motor processing in primary visual cortex"

startover

%% Or load dataset:
load Dataset2_1.mat



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
cell_IDs{end+1}     =     '20102018081611378'; %new paper example

 
% %Other examples:
% cell_IDs{end+1}     = '20122018081611121'; %auditory unrelated to movement
% cell_IDs{end+1}     = '20102018081611325';
% cell_IDs{end+1}     = '10122019041031329'; %Auditory no motor example cell glm
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
params.exportfig    = 1;
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




