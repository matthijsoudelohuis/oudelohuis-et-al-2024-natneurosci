% Oude Lohuis et al. 2023
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

%% Load opto only:
params.animals                  = {'1008' '1009' '1011'};
params.experiment               = 'OptoOnly';
params.AlignOn                  = 'photostimStart';

%% Load data:
[Data]          = MOL_GetData('E:','CHDET',params.experiment,params.animals,[],{'sessionData' 'trialData' 'lfpData'});

sessionData     = Data.sessionData;
trialData       = Data.trialData;
lfpData         = Data.lfpData;

%% Filter out V1 channels:
keepidx = strcmp(lfpData.area,params.area);
lfpfields = fieldnames(lfpData);
for iF = 1:length(lfpfields)
    lfpData.(lfpfields{iF}) = lfpData.(lfpfields{iF})(keepidx,:);
end
fprintf('Filtered %d V1 channels\n',length(lfpData.session_ID));

%% Flip signs:
lfpData.sortedChannelDepth = -lfpData.sortedChannelDepth;
if params.flipsignLFP
    lfpData.signal = cellfun(@(x) -x, lfpData.signal, 'UniformOutput',false);
end

%% Save dataset: % large file with all the lfp samples
save('DatasetS9_3.mat','params','sessionData','trialData','lfpData')

%% Or start script from saved dataset:
load DatasetS9_3.mat




%% Init output: interpolated CSD matrix:
nSessions                       = length(sessionData.session_ID);
outMeanCSD                      = NaN(length(params.finalYaxis),nTimebins,nSessions);
outMeanERP                      = NaN(length(params.finalYaxis),nTimebins,nSessions);

%% For each session compute CSD:
for iSes = 1:nSessions
    fprintf('Computing CSD for session %d/%d\n',iSes,nSessions)
    
    %Get the relevant data for each session individually with the actual lfp data:
    [tempsessionData,temptrialData,templfpData]        = MOL_getTempPerSes(sessionData.session_ID(iSes),sessionData,trialData,lfpData);

    if ~isempty(templfpData.signal) %&& any(templfpData.lfpgood)
        
        switch tempsessionData.(sprintf('Probe%d_Config',unique(templfpData.probeIdx))){1}
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
        
        params.Ylim_all                    	= [min(templfpData.sortedChannelDepth)-25 max(templfpData.sortedChannelDepth)+25];
        params.nChannels                    = length(templfpData.ch);
        
        templfpData.ts                      = repmat(templfpData.t_start(1):1/templfpData.fs(1)*1e6:templfpData.t_end(1),1,1);
        
        templfpData.sortedChannelNum        = mod(templfpData.sortedChannelNum-1,params.nChannels)+1;
        templfpData.sortedChannelidx        = mod(templfpData.sortedChannelidx-1,params.nChannels)+1;

        templfpData.sortedSignal            = cell2mat(templfpData.signal);
        templfpData.sortedSignal            = templfpData.sortedSignal(templfpData.sortedChannelidx,:);
        
        params.ChannelSel                   = 1:params.nChannels;
        
        events_ts                           = temptrialData.(params.AlignOn);
        
        [meancsd,meanerp]                   = MOL_CSD(params, events_ts, templfpData);
        
%         idx_ch                              = templfpData.sortedChannelDepth>params.chdepthmin-50        & templfpData.sortedChannelDepth<params.chdepthmax+50;
        idx_ch                              = templfpData.sortedChannelDepth>-params.chdepthmax-50       & templfpData.sortedChannelDepth<-params.chdepthmin+50;
        
        idx_out                             = params.finalYaxis>min(templfpData.sortedChannelDepth)-50     & params.finalYaxis<max(templfpData.sortedChannelDepth)+50;
        if params.interpolate
            outMeanCSD(idx_out,:,iSes)                        = tointerpol2(meancsd(idx_ch,:),templfpData.sortedChannelDepth(idx_ch),params.finalYaxis(idx_out));
            outMeanERP(idx_out,:,iSes)                        = tointerpol2(meanerp(idx_ch,:),templfpData.sortedChannelDepth(idx_ch),params.finalYaxis(idx_out));
        end
    end
end


%% CSD at different depths (not in MS)
sesfilter = strcmp(sessionData.PhotostimArea,'A1') & sessionData.Photostimpower>=1 & sessionData.PhotostimFreq==20;
fprintf('Average CSD Fiber at A1: n=%d sessions\n',sum(sesfilter))

%compute mean only to set the colorbarscale:
meancsd                             = squeeze(nanmean(outMeanCSD(:,:,sesfilter),3));
meanerp                             = squeeze(nanmean(outMeanERP(:,:,sesfilter),3));
params.cscale                       = [-max(max(max(meancsd)))*0.92 max(max(max(meancsd)))*0.92];

csdfig                              = figure; set(gcf,'units','normalized','Position',[0.6 0.4 0.25 0.4],'color','w');

imagesc(params.xtime,params.finalYaxis,meancsd,params.cscale); hold on;
xlabel('Time from stimulus (s)','FontSize', 15)
% set(gca,'YDir','reverse') %flip the direction of imagesc
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
%To plot the ERP over it, compute an offset:
offsetmat = repmat(params.finalYaxis',1,length(params.xtime));
plot(params.xtime,meanerp*50e4 + offsetmat,'k','LineWidth',0.5); hold on; %plot mean erp with offset
if strcmp(params.experiment,'OptoOnly')
    plot([0 10e3],[-5 -5],'b','LineWidth',7)
    plot([50e3 60e3],[-5 -5],'b','LineWidth',7)
end
%Figure make up:
title('Fiber over A1','FontSize',15)
%Figure make up:
ylabel('Depth','FontSize', 15)
plot([0 0],ylim,'k','LineWidth',1);
xlim([params.xtime(1) params.xtime(end)]);
set(gca,'XTick',params.xticks,'XTickLabels',params.xticklabels,'FontSize',15)
set(gca,'YTick',params.yticks,'YTickLabels',params.yticklabels,'FontSize',15)
xlabel('Time from laser onset (ms)','FontSize', 15)
ylim([-params.chdepthmax -params.chdepthmin])

%% Fig S9k Compute mean ERP response as a function of photostimulation power: 
sesfilter                       = strcmp(sessionData.PhotostimArea,'A1') & sessionData.PhotostimFreq==20 & strcmp(sessionData.Rec_datetime,'2019-03-06_12-04-15');

nSessions                       = sum(sesfilter);
powers                          = sessionData.Photostimpower(sesfilter);

idx_time                        = params.xtime>5e3 & params.xtime<20e3;

params.colors_power             = getPyPlot_cMap('gnuplot2',10);

meanerp_resp                    = squeeze(nanmean(outMeanERP(:,idx_time,sesfilter),2));

figure;  set(gcf,'units','normalized','Position',[0.6 0.4 0.14 0.2],'color','w'); hold all;
mean_erp_amp = nanmean(abs(meanerp_resp),1);
plot(powers',mean_erp_amp,'-','Color','k','LineWidth',2)
for i = 1:length(powers)
    plot(powers(i),mean_erp_amp(i),'.','Color',params.colors_power(i,:),'MarkerSize',30)
end
ylim([0 0.5e-4])
xlim([-0.15 10])
ylabel('|ERP amplitude| (\muV)')
xlabel('Fiber power (mW)')
set(gca,'XTick',powers,'YTick',[0 25 50]*1e-6,'YTickLabels',[0 25 50])
set(gca,'XTick',[0 1 5 10],'YTick',[0 25 50]*1e-6,'YTickLabels',[0 25 50])
MOL_prepfigAI

fprintf('Average CSD Fiber at A1: n=%d sessions\n',nSessions)
export_fig(fullfile(params.savedir,'ERP_laserpower'),'-eps');


