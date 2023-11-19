%% Oude Lohuis et al. 2023 Nat Neurosci
% Triple dissociation of auditory, visual, and motor processing in primary visual cortex
% MOL (C) 2023

% Script to plot raw trace after muscimol injection

%Data not included in the repository!

%% 
params              = params_histresponse;
params.conv_twin    = 10e6;        %Window size for smoothing
params.conv_sigma   = 2e6;        %sd of gaussian window for smoothing

%% Get arguments:
params.RootDir      = 'F:';
params.LocalDir     = 'E:';

RootDir             = 'F:\Data\CHDET\RawData\2025\2020-01-23_10-36-27\';

%Save figure as given by user input
params.savedir       = 'E:\Documents\PhD\Figures\Project CHDET\Results - auV1\7MuscimolA1';

Alpha               = 0.3;
params.Alpha        = 0.3;

params.examplech    = 60;
params.examplech    = 3;

params.nthresholds = 5;

%% Get the data: %do not do any filtering or resampling, get raw data
%Parameters for Butterworth filter
params.UseButter               = 0;
%Parameters for Kaiser filter
params.UseKaiser               = 0;
%Parameters for resampling
params.UseResampling           = 0;    %Resample data to lower frequency

%% Get the events file:
filename                    = fullfile(RootDir,'CheckersMuscimol','events.nev');
sessionData                 = struct();
sessionData.Experiment      = 'CheckersMuscimol';
[sessionData, trialData]    = MOL_Preprocessing_Checkerboard(sessionData,RootDir);

trialData.stimStart         = trialData.stimChange;
sessionData.MuscimolStart = sessionData.MuscimolStart-10e6;
TimeWindow                  = [sessionData.t_start sessionData.t_stop];

%% Get the data:
lfpData                         = MOL_extractLFP(fullfile(RootDir,'RawData'),params.examplech,TimeWindow,params);

%Apply high pass filter:
params.hp_butter                = 500;  %Low pass filter (Hz) (Butter)
params.ord_butter               = 4;    %Butterworth filter order
[B_butt,A_butt]                 = butter(params.ord_butter,params.hp_butter/(lfpData.fs(1)/2),'high');

lfpData.hpsignal{1}             = filtfilt(B_butt,A_butt,lfpData.signal{1});

lfpData.ts = lfpData.t_start:(1e6/lfpData.fs):lfpData.t_end;

exampletrials = [];

%% Filter out artefacts:
baseline            = lfpData.ts<sessionData.MuscimolStart;
idx                 = find(abs(lfpData.hpsignal{1})>13*std(lfpData.hpsignal{1}(baseline)));

h = repmat(idx,11,1);
j = repmat(transpose(-5:5),1,length(idx));
g = h + j;
g = unique(g);
lfpData.hpsignal{1}(g) = 0;

%% Plot the high pass filtered trace:
figure; 
set(gcf,'units','normalized','Position',[0.2 0.25 0.4 0.55],'color','w')
selectedpiece           = find(lfpData.ts>TimeWindow(1) & lfpData.ts<TimeWindow(end));

% subsample:
selectedpiece = selectedpiece(1:50:end);
plot(lfpData.ts(selectedpiece),lfpData.hpsignal{1}(selectedpiece),'k'); hold all;
ymax    = max(lfpData.hpsignal{1}(selectedpiece))*1.1;
ymin    = min(lfpData.hpsignal{1}(selectedpiece))*1.1;

plot([sessionData.MuscimolStart sessionData.MuscimolStart],[ymin*1.2 ymax*1.2],'k','LineWidth',2); hold all;

ylim([ymin*1.2 ymax*1.2]);
xlim(TimeWindow);

interval = 5;
nintervals = 5;
xticks = [sessionData.MuscimolStart:interval*60e6:sessionData.MuscimolStart+nintervals*60e6*interval];
xticklabels = 0:interval:nintervals*interval;

set(gca,'linewidth',3)
set(gca,'XTick',xticks,'XTickLabel',xticklabels,'FontSize',25)
xlabel('Time after injection (min)','FontSize',30)
ylabel('Voltage (uv)','FontSize',30)

export_fig(fullfile(params.savedir,strcat('Ch7_Raw_Trace','.eps')),'-eps');

%% plot MUA over time:
baseline            = lfpData.ts<sessionData.MuscimolStart;
idx                 = find(abs(lfpData.hpsignal{1})>params.nthresholds*std(lfpData.hpsignal{1}(baseline)));
idx                 = idx(diff(idx)~=1);

edges               = lfpData.t_start:params.binsize:lfpData.t_end;
hist_mat            = histc(lfpData.ts(idx),edges) * 1e6/params.binsize;

N                   = params.conv_twin/params.binsize;
alpha               = ((N-1)/(params.conv_sigma/params.binsize))/2; %? = (N – 1)/(2?)
win                 = gausswin(N,alpha); %convolution with gaussian
win                 = win/sum(win); %normalized

%Smooth either the total or the individual trials:
hist_mat            = padarray(hist_mat,[0 round(length(win)/2)],'symmetric','both'); %pad the array on both sides for convolution
hist_mat            = conv(hist_mat,win,'valid'); %Take only the valid overlapping center of convolution
hist_mat            = hist_mat(1:length(edges)); %slight correction to get same size (edge vs convolution)

hist_mat = hist_mat/(hist_mat(1));
figure; hold all; set(gcf,'units','normalized','Position',[0.2 0.2 0.22 0.11],'color','w');

%Make figure:
selectedpiece           = edges>TimeWindow(1) & edges<TimeWindow(end);

plot(edges,hist_mat,'k','LineWidth',1); hold all;
ymax    = max(hist_mat(selectedpiece))*1.1;
ymin    = min(hist_mat(selectedpiece))*1.1;

ylim([ymin*1.2 ymax*1.2]);
xlim(TimeWindow);

interval = 5;
nintervals = 5;
xticks = [sessionData.MuscimolStart:interval*60e6:sessionData.MuscimolStart+nintervals*60e6*interval];
xticklabels = 0:interval:nintervals*interval;

set(gca,'XTick',xticks,'XTickLabel',xticklabels,'YTick',[1 2])

ylim([0 3])
ylabel('Norm. AC MUA')
xlabel('Time in session')

export_fig(fullfile(params.savedir,'Rateses_AC'),'-eps','-nocrop')
