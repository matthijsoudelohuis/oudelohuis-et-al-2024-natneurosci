%% Oude Lohuis et al. 2023 Nat Neurosci
% Triple dissociation of auditory, visual, and motor processing in primary visual cortex
% MOL (C) 2023

% This script describes shows the motion energy overlaid with the video
% as reported in Oude Lohuis et al. 2023
% Extended data figure 2

%# DATA NOT INCLUDED IN THE REPOSITORY! 

%% Parameters:
params.shiftNframes     = 2;
params.spat_filter      = 1;
params.savedir          = 'E:\Documents\PhD\Figures\Project CHDET\Results - auV1\13MovieBehavior\VideoProjection';

%% Example NE mouse:
params.videoName        = 'VT1.0001';
params.rootDir          = 'G:';
params.animal           = '1012';
params.session          = '2019-04-11_11-12-47';
params.experiment       = 'ChangeDetectionConflictDecor';
params.shiftNframes     = 2;
% params.t_exframes       = [-400e3 80e3 200e3 1400e3];
params.t_exframes       = [0e3 80e3 200e3 520e3 1400e3];

params.clims            = [150 650];
params.framenumber      = 1042;

plotMotionOverlay(params);

%% Active MST example:
params.videoName        = 'VT1.0001';
params.rootDir          = 'F:';
params.animal           = '2031';
params.session          = '2020-01-24_15-19-14';
params.experiment       = 'ChangeDetectionConflict';
params.shiftNframes     = -2;
% params.t_exframes       = [-400e3 80e3 200e3 1400e3];
params.t_exframes       = [0e3 80e3 200e3 520e3 1400e3];

params.clims            = [130 500];
params.framenumber      = 1807;

plotMotionOverlay(params);

%% Passive MST example:
params.videoName        = 'VT1.0001';
params.rootDir          = 'F:';
params.animal           = '2045';
params.session          = '2021-04-29_15-55-10';
params.experiment       = 'ChangeDetectionConflict';
params.shiftNframes     = 2;
% params.t_exframes       = [-400e3 80e3 200e3 1400e3];
params.t_exframes       = [0e3 120e3 320e3 520e3 1400e3];

params.clims            = [130 350];
params.framenumber      = 1042;

plotMotionOverlay(params);

