% Script to plot modulation by CamKIIa neurons with opto:
% Oude Lohuis et al. 2023

% Data not included in the repository!!! (too large)

startover

%% Parameters:
params.linewidth        = 0.25; %stroke weight

exampletrials           = [2 5]; %select 4 trials out 10
params.alpha            = 1; %for print
pulselength             = 10e3; %individual pulse length
fs                      = 32000; %sampling frequency of neuralynx
params.savedir          = 'E:\Documents\PhD\Figures\Project CHDET\Results - auV1\4CaMKIIa_excitation\';% Output dir settings:

%% Session and channel:
RootDir                 = 'G:\';
exsession               = {'2019-03-06_12-04-21'};
selectedchannel         = 44;
params.UseResampling    = 0;    %Resample data to lower frequency

%% Get session and trial data:
[Data] = MOL_GetData('E:','CHDET','OptoOnly',[],exsession,{'sessionData' 'trialData'});
sessionData             = Data.sessionData;
trialData               = Data.trialData;
RawDataDir              = fullfile(RootDir,'Data','CHDET','RawData',sessionData.mousename{1},sessionData.Rec_datetime{1},'RawData');

%% Get the lfp data:
TimeWindow              = [sessionData.t_start sessionData.t_stop];
lfpData                 = MOL_extractLFP(RawDataDir,selectedchannel,TimeWindow,params);
lfpData.ts              = linspace(sessionData.t_start,sessionData.t_stop,length(lfpData.signal{1}));

%alternatively: load lfpData from savefile:
load('DatasetX_lfpData.mat','lfpData');

%% Plot for a few trials:
figure; 
example_tstart          = trialData.photostimStart(exampletrials(1))-1e6;
example_tstop           = trialData.photostimStart(exampletrials(2))+2e6;
set(gcf,'units','normalized','Position',[0.2 0.5 0.37 0.257],'color','w')
selectedpiece           = lfpData.ts>example_tstart & lfpData.ts<example_tstop;
plot(lfpData.ts(selectedpiece),lfpData.signal{1}(selectedpiece),'k','LineWidth',params.linewidth); hold all;
ymax    = max(lfpData.signal{1}(selectedpiece))*1.1;
ymin    = min(lfpData.signal{1}(selectedpiece))*1.1;

for iT = exampletrials(1):exampletrials(2)
    optostart   = trialData.photostimStart(iT);
    optoend     = trialData.photostimEnd(iT);
    X = [optostart optoend optoend optostart];
    Y = [ymax-0.0001 ymax-0.0001 ymax ymax];
    patch(X, Y,'blue','EdgeColor','blue','FaceAlpha',params.alpha,'EdgeAlpha',params.alpha); hold on;
end

set(gca,'visible','off') %remove axes, only trace

ylim([ymin*1.2 ymax*1.2]);
xlim([example_tstart example_tstop]);

sc = scalebar;
sc.XLen = 1e6;
sc.XUnit = 'sec';
sc.YLen = 0.5e-3;
sc.YUnit = 'V';

print(gcf,fullfile(params.savedir,strcat('Ch',num2str(selectedchannel),'_Raw_Trace_','.pdf')),'-dpdf','-bestfit');
export_fig(fullfile(params.savedir,strcat('Ch',num2str(selectedchannel),'_Raw_Trace_','.png')),'-png');

%% Make figure of close up:
figure;
set(gcf,'units','normalized','Position',[0.2 0.1 0.5 0.7],'color','w')

preposttime         = 0.5e6;
exampletrial        = 2;
selectedpiece       = lfpData.ts>trialData.photostimStart(exampletrial)-preposttime & lfpData.ts<trialData.photostimEnd(exampletrial)+preposttime;

plot(lfpData.ts(selectedpiece),lfpData.signal{1}(selectedpiece),'k','LineWidth',params.linewidth); hold all;
ymax = max(lfpData.signal{1}(selectedpiece))*1.1;
ymin = min(lfpData.signal{1}(selectedpiece))*1.1;

for i = 1:sessionData.PhotostimFreq
    pulsestart = trialData.photostimStart(exampletrial)+1e6/sessionData.PhotostimFreq*(i-1);
    X = [pulsestart pulsestart+pulselength pulsestart+pulselength pulsestart];
%     Y = [ymax-0.001 ymax-0.001 ymax ymax];
    Y = [ymax-0.0001 ymax-0.0001 ymax ymax];
    patch(X, Y,'blue','EdgeColor','blue','FaceAlpha',params.alpha,'EdgeAlpha',params.alpha); hold on;
end

set(gca,'visible','off')

ylim([ymin*1.1 ymax*1.1]);
set(gca,'linewidth',3)
print(gcf,fullfile(params.savedir,strcat('Ch',num2str(selectedchannel),'_Pulse_Closeup_',num2str(sessionData.Photostimpower),'mW','.pdf')),'-dpdf','-bestfit');
export_fig(fullfile(params.savedir,strcat('Ch',num2str(selectedchannel),'_Pulse_Closeup_',num2str(sessionData.Photostimpower),'mW','.png')),'-png');

%% plot MUA over time:
nthresholds = 4;

idx                 = find(abs(lfpData.signal{1})>nthresholds*std(lfpData.signal{1}));
idx                 = idx(diff(idx)~=1);

params              = params_histresponse_coding(params);
params.conv_sigma           = 0.1e6;        %sd of gaussian window for smoothing

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

figure;
set(gcf,'units','normalized','Position',[0.2 0.1 0.5 0.7],'color','w')
selectedpiece           = edges>example_tstart & edges<example_tstop;
% area(edges(selectedpiece),hist_mat(selectedpiece)); hold all;

area(edges(selectedpiece),hist_mat(selectedpiece),'FaceColor','k'); hold all;
% plot(edges(selectedpiece),hist_mat(selectedpiece),'k'); hold all;
ymax    = max(hist_mat(selectedpiece))*1.1;
ymin    = min(hist_mat(selectedpiece))*1.1;

for iT = exampletrials(1):exampletrials(2)
    optostart   = trialData.photostimStart(iT);
    optoend     = trialData.photostimEnd(iT);
    X = [optostart optoend optoend optostart];
    Y = [ymax-0.5 ymax-0.5 ymax ymax];
    patch(X, Y,'blue','EdgeColor','blue','FaceAlpha',params.alpha,'EdgeAlpha',params.alpha); hold on;
end

ylim([ymin*1.2 ymax*1.2]);
xlim([example_tstart example_tstop]);
set(gca,'linewidth',3)
print(gcf,fullfile(params.savedir,strcat('Ch',num2str(selectedchannel),'_MUA_Trace','.pdf')),'-dpdf','-bestfit');
export_fig(fullfile(params.savedir,strcat('Ch',num2str(selectedchannel),'_MUA_Trace','.png')),'-png');

