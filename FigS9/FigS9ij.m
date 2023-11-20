%% Oude Lohuis et al. 2024 Nat Neurosci
startover

%% Parameters:

%General parameters
params.t_pre                    = -0.02e6; %All time in microseconds
params.t_post                   = .1e6;  %All time in microseconds

% CSD parameters:
params.flipsignLFP              = 1;    %Whether to flip the LFP because NLX has inversed polarity (so flip again to make comparable to others)
params.csdmethod                = 'NicholsonFreeman'; %'NicholsonFreeman'; %{'NicholsonFreeman' 'DoubleDiff'}
params.colormap                 = 'parula'; %redblue or parula
params.sinkhot                  = 1;    %Whether sink is shown as red hot colors (2 conventions exist in csd world)
params.conductivity             = 0.4;  %in Siemens m-1
% There are filtering before making the CSD:
% Parameters for Butterworth filter
params.UseButter                = 0;
params.lp_butter                = 300;  %Low pass filter (Hz) (Butter)
params.ord_butter               = 4;   %Butterworth filter order

%Parameters for Kaiser filter
params.UseKaiser                = 0;
params.lp_kaiser                = 80;  %Low pass filter (Hz) (Kaiser)
params.hp_kaiser                = 0.1;   %High pass filter (Hz) (Kaiser)
params.dp_kaiser                = 1;   %Transition bandwidth (Kaiser only)

params.area                     = {'V1'};

% params.chdepthmin               = 25;
% params.chdepthmax               = -1050;

params.chdepthmin               = -25;
params.chdepthmax               = 1000;

params.interpolate              = 1;
params.resolution               = 20;
params.finalYaxis               = params.chdepthmin:params.resolution:params.chdepthmax;
params.finalYaxis               = -params.chdepthmax:params.resolution:-params.chdepthmin;

params.yticks                   = [-1000:200:-0];
params.yticklabels              = -params.yticks;

params.xticks                   = params.t_pre:10e3:params.t_post;
params.xticklabels              = params.xticks*1e-3;

params.fs                       = 1024;

params.xtime                    = params.t_pre:1e6/params.fs:params.t_post(1) - 1e6/params.fs;
nTimebins                       = length(params.xtime);

params.savedir                  = 'E:\Documents\PhD\Figures\Project CHDET\Results - auV1\4CaMKIIa_excitation\';% Output dir settings:

%% Run this block for the optical photostim at A1 at 10 mW and show CSD in V1 
sestitle                        = 'Fiber over A1';
params.experiment               = 'OptoOnly';
params.exsession                = '2019-03-06_12-04-22';
params.exsession                = '2019-03-06_12-04-21';
params.AlignOn                  = 'photostimStart';
movartefac                      = [21 20 12]; %Channels with movement artefact (interpolate out of precaution)

% Load data:
[Data]          = MOL_GetData('E:','CHDET',params.experiment,[],params.exsession,{'sessionData' 'trialData' 'lfpData'});
% [Data]          = MOL_GetData('E:','CHDET',{'OptoOnly'},{'1009'},{'2019-03-07_12-46-61'},{'sessionData' 'trialData' 'lfpData'});

sessionData     = Data.sessionData;
trialData       = Data.trialData;
lfpData         = Data.lfpData;

keepidx = strcmp(lfpData.area,params.area);
lfpfields = fieldnames(lfpData);
for iF = 1:length(lfpfields)
    lfpData.(lfpfields{iF}) = lfpData.(lfpfields{iF})(keepidx,:);
end

lfpData.sortedChannelDepth = -lfpData.sortedChannelDepth;

if params.flipsignLFP
    lfpData.signal = cellfun(@(x) -x, lfpData.signal, 'UniformOutput',false);
end

%% Save dataset:
save('DatasetS9_1.mat','params','sessionData','trialData','lfpData')

%% Or start script from saved dataset:
load DatasetS9_1.mat


%% Run this block for the comparative visual stimulation:
sestitle                        = 'Checkers';
params.experiment               = 'CheckerboardReversals';
params.exsession                = '2019-03-06_12-04-15';
params.AlignOn                  = 'stimStart';
movartefac                      = [21 20 12]; %Channels with movement artefact (interpolate out of precaution)

%% Load data:
[Data]          = MOL_GetData('E:','CHDET',params.experiment,[],params.exsession,{'sessionData' 'trialData' 'lfpData'});
% [Data]          = MOL_GetData('E:','CHDET',{'OptoOnly'},{'1009'},{'2019-03-07_12-46-61'},{'sessionData' 'trialData' 'lfpData'});

sessionData     = Data.sessionData;
trialData       = Data.trialData;
lfpData         = Data.lfpData;

keepidx = strcmp(lfpData.area,params.area);
lfpfields = fieldnames(lfpData);
for iF = 1:length(lfpfields)
    lfpData.(lfpfields{iF}) = lfpData.(lfpfields{iF})(keepidx,:);
end

lfpData.sortedChannelDepth = -lfpData.sortedChannelDepth;

if params.flipsignLFP
    lfpData.signal = cellfun(@(x) -x, lfpData.signal, 'UniformOutput',false);
end

%% Save dataset:
save('DatasetS9_2.mat','params','sessionData','trialData','lfpData')

%% Or start script from saved dataset:
load DatasetS9_2.mat


%% Make CSD and ERP of the data:
lfpData.lfpgood(lfpData.hasMovArtefacts==1)         = 0;

switch sessionData.(sprintf('Probe%d_Config',unique(lfpData.probeIdx))){1}
    case 'A4x8-5mm-100-200-177-CM32'
        params.nChannels           = 32;
        params.nShanks             = 4;
        params.nChannelsShank      = 8;
        params.intersitedistance   = 100;
    case 'A1x32-Poly2-5mm-50s-177-CM32'
        params.nChannels           = 32;
        params.nShanks             = 2;
        params.nChannelsShank      = 16;
        params.intersitedistance   = 50;
    case 'A1x64-Poly2-6mm-23s-160'
        params.nChannels           = 64;
        params.nShanks             = 2;
        params.nChannelsShank      = 32;
        params.intersitedistance   = 46;
end

lfpData.ts                      = repmat(lfpData.t_start(1):1/lfpData.fs(1)*1e6:lfpData.t_end(1),1,1);

lfpData.sortedChannelNum        = mod(lfpData.sortedChannelNum-1,params.nChannels)+1;
lfpData.sortedChannelidx        = mod(lfpData.sortedChannelidx-1,params.nChannels)+1;

lfpData.sortedSignal            = cell2mat(lfpData.signal);
lfpData.sortedSignal            = lfpData.sortedSignal(lfpData.sortedChannelidx,:);

% lfpData.sortedisgood(64-15) = 0;
% lfpData.sortedisgood(43) = 1;
lfpData.sortedisgood(movartefac) = 0;
% lfpData.sortedisgood = 0;

params.ChannelSel                   = 1:params.nChannels;

events_ts                           = trialData.(params.AlignOn);
% events_ts                           = trialData.photostimStart(2:end);

[meancsd,meanerp]                   = MOL_CSD(params, events_ts, lfpData);

% Figures and analyses:
%initialize figure
csdfig                              = figure; set(gcf,'units','normalized','Position',[0.6 0.4 0.25 0.4],'color','w');

params.cscale                       = [-max(max(max(meancsd)))*0.92 max(max(max(meancsd)))*0.92];

% params.cscale                       = [-1 1];

%Plot the CSD:
% imagesc(params.xtime,params.finalYaxis,meancsd,params.cscale); hold on;
imagesc(params.xtime,lfpData.sortedChannelDepth,meancsd,params.cscale); hold on;
xlabel('Time from stimulus (s)','FontSize', 15)
set(gca,'YDir','normal') %flip the direction of imagesc

switch params.colormap
    case 'redblue'
        h       = coolwarm + repmat(0.1-abs(linspace(-0.1,.1,64))',1,3);
        h(h>1)  = 1;
        if params.sinkhot
            h = flipud(h);
        end
        colormap(h);
    case 'parula'
        if params.sinkhot
            colormap(flipud(parula));
        else                     colormap(parula);
        end
end

colorbar();
%To plot the ERP over it, compute an offset:
offsetmat = repmat(lfpData.sortedChannelDepth,1,length(params.xtime));
plot(params.xtime,meanerp*50e4 + offsetmat,'k','LineWidth',0.5); hold on; %plot mean erp with offset
% plot(params.xtime,meanerp(movartefac,:)*50e4 + offsetmat(movartefac,:),'k','LineWidth',5); hold on; %plot mean erp with offset
if strcmp(params.experiment,'OptoOnly')
    plot([0 10e3],[-5 -5],'b','LineWidth',7)
    plot([50e3 60e3],[-5 -5],'b','LineWidth',7)
end
%Figure make up:
title(sestitle,'FontSize',15)
%Figure make up:
ylabel('Depth','FontSize', 15)
plot([0 0],ylim,'k','LineWidth',1);
xlim([params.xtime(1) params.xtime(end)]);
set(gca,'XTick',params.xticks,'XTickLabels',params.xticklabels,'FontSize',15)
set(gca,'YTick',params.yticks,'YTickLabels',params.yticklabels,'FontSize',15)
xlabel('Time from laser onset (ms)','FontSize', 15)
ylim([-params.chdepthmax -params.chdepthmin])

