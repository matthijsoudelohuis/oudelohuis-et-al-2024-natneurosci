%% Parameter settings:
params                      = params_histresponse(); % All time is in microseconds
params.Experiments          = {'ChangeDetectionConflictDecor' 'VisOnlyTwolevels' 'ChangeDetectionConflict' }; %Which versions of the task to load data from
params.ExperimentLabels     = {'NE' 'UST' 'MST'}; %Which versions of the task to load data from

% params.exSession            = '2019-03-07_12-46-53';
% params.exSession            = '2019-04-10_11-43-25';
% params.exSession            = '2019-10-22_17-57-16';
params.exSession            = '2019-10-24_18-21-04'; %Supp figure example
% params.exSession            = '2019-07-16_15-50-06';

params.savedir              = 'E:\Documents\PhD\Figures\Project CHDET\Results - auV1\15LayerEstimation';

params.colormap_CSD         = 'parula';
params.sinkhot              = 1;

params.depthcorrection      = -50;

params.t_pre                = -0.1e6; %All time in microseconds
params.t_post               = .4e6;  %All time in microseconds

params.xticks               = params.t_pre:1e5:params.t_post;
params.xticklabels          = params.xticks*1e-3;

params.fs                   = 32e3;

params.xtime                = params.t_pre:1e6/params.fs:params.t_post(1) - 1e6/params.fs;
params.nTimebins            = numel(params.xtime);

%% Get data:
[Data] = MOL_GetData('E:','CHDET',params.Experiments,[],params.exSession,{'sessionData' 'trialData' 'lfpData_lazy'});
sessionData     = Data.sessionData;
trialData       = Data.trialData;
lfpData         = Data.lfpData;

%% Filter channels on area:
idx             = ismember(lfpData.area,{'V1'});
lfpFields     = fieldnames(lfpData);
for iF = 1:length(lfpFields)
    lfpData.(lfpFields{iF}) = lfpData.(lfpFields{iF})(idx,:);
end

%%  MUA power:
figure; hold all; set(gcf,'units','normalized','Position',[0.12 0.45 0.41 0.34],'color','w');

subplot(1,2,1);

%Get power estimate for different frequency bands:
[tempsessionData,temptrialData,templfpData] = MOL_getTempPerSes(sessionData.session_ID(1),sessionData,trialData,lfpData);
templfpData.sortedChannelDepth = templfpData.sortedChannelDepth + params.depthcorrection;

templfpData.sortedChannelDepth = -templfpData.sortedChannelDepth;

% templfpData.HF_PWR(templfpData.sortedisgood==1)
plot(templfpData.HF_PWR(templfpData.sortedisgood==1),templfpData.sortedChannelDepth(templfpData.sortedisgood==1),'k.-','LineWidth',2,'MarkerSize',25)
xlabel('Normalized MUA power')
ylim([-1000 0])
set(gca,'YDir','normal')

%Panel 3: CSD power
subplot(1,2,2)
meancsd             = templfpData.CheckerCSD;
meanerp            = templfpData.CheckerERP;

params.cscale               = [-max(max(meancsd))*0.95 max(max(meancsd))*0.95];

imagesc(params.xtime,templfpData.sortedChannelDepth,meancsd,params.cscale); hold on;
switch params.colormap_CSD
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

offsetmat = repmat(templfpData.sortedChannelDepth,1,length(params.xtime));
% plot(params.xtime,meanerp*20e4 + offsetmat,'k','LineWidth',0.5); hold on;
idx_subsample = 1:10:params.nTimebins;
plot(params.xtime(idx_subsample),meanerp(:,idx_subsample)*10e4 + offsetmat(:,idx_subsample),'k','LineWidth',0.5); hold on;

ylabel('Channel','FontSize', 15)
plot([0 0],ylim,'k','LineWidth',1);
xlim([params.xtime(1) params.xtime(end)]);
set(gca,'XTick',params.xticks,'XTickLabels',params.xticklabels)
xlabel('Time from stimulus (ms)','FontSize', 15)

ylim([-1000 0])
set(gca,'YDir','normal')

% export_fig(fullfile(params.savedir,sprintf('ExSession_%s',tempsessionData.session_ID{1})),'-eps','-nocrop')
export_fig(fullfile(params.savedir,sprintf('CSD_MUA_ExSession_%s',tempsessionData.session_ID{1})),'-eps')

%% SESSION AVERAGES: 



%% Get data for example session:
% [Data] = MOL_GetData('E:','CHDET',params.Experiments,{'2003' '2011' '2012' '2030' '2031'},[],{'sessionData' 'trialData' 'lfpData_lazy'});
[Data] = MOL_GetData('E:','CHDET',params.Experiments,[],[],{'sessionData' 'trialData' 'lfpData_lazy'});
% [Data] = MOL_GetData('E:','CHDET',params.Experiments,{},[],{'sessionData' 'trialData' 'lfpData_lazy'});
sessionData     = Data.sessionData;
trialData       = Data.trialData;
lfpData         = Data.lfpData;

%% Filter channels on area:
idx             = ismember(lfpData.area,{'V1'});
lfpFields     = fieldnames(lfpData);
for iF = 1:length(lfpFields)
    lfpData.(lfpFields{iF}) = lfpData.(lfpFields{iF})(idx,:);
end

%% Get only sessions in V1 with L5 max and CSD data:
sesids                          = unique(lfpData.session_ID(lfpData.L5_maxMUA==1 & strcmp(lfpData.area,'V1')));

[sessionData,trialData,lfpData] = MOL_getTempPerSes(sesids,sessionData,trialData,lfpData);
fprintf('Dataset: %d sessions, %d trials, %d channels\n',length(sessionData.session_ID),length(trialData.session_ID),length(lfpData.session_ID));

nSessions                       = length(sessionData.session_ID);

%% Copy hf power to sorted field:
for iSes = 1:nSessions
    idx_ses                         = strcmp(lfpData.session_ID,sessionData.session_ID(iSes));
    if all(isnan(lfpData.sortedHF_PWR(idx_ses)))
        lfpData.sortedHF_PWR(idx_ses) = lfpData.HF_PWR(idx_ses);
    end
end

%% Apply depth correction:
lfpData.sortedChannelDepth      = lfpData.sortedChannelDepth + params.depthcorrection;
lfpData.ChannelY                = lfpData.ChannelY + params.depthcorrection;

%and invert depth:
lfpData.sortedChannelDepth      = -lfpData.sortedChannelDepth;
lfpData.ChannelY                = -lfpData.ChannelY;

%%

figure; hold all; set(gcf,'units','normalized','Position',[0.02 0.45 0.16 0.3],'color','w');

params.colormap_probes = getPyPlot_cMap('brg',nSessions);

for iSes    = 1:nSessions
    y_min   = min(lfpData.ChannelY(strcmp(lfpData.session_ID,sessionData.session_ID(iSes))));
    y_max   = max(lfpData.ChannelY(strcmp(lfpData.session_ID,sessionData.session_ID(iSes))));
    plot([iSes iSes],[y_min y_max],'-','LineWidth',0.6,'Color',params.colormap_probes(iSes,:))
    
end
xlabel('Session')
plot([0 nSessions],[-1000 -1000],'k:','LineWidth',0.5)
plot([0 nSessions],[0 0],'k:','LineWidth',0.5)
plot([0 nSessions],[1000 1000],'k:','LineWidth',0.5)
xlim([0.5 nSessions+0.5])
ylim([-100 1100])
set(gca,'YDir','reverse','XTick',[])

ylim([-1100 100])
set(gca,'YDir','normal','XTick',[])

%% Init output: interpolated CSD matrix:

params.interpolate              = 1;
% params.chdepthmin               = 0;
% params.chdepthmax               = 1000;

params.chdepthmin               = -1000;
params.chdepthmax               = 0;

params.resolution               = 20;
params.finalYaxis               = params.chdepthmin:params.resolution:params.chdepthmax;

nSessions                       = length(sessionData.session_ID);
nTimebins                       = length(params.xtime);

outMeanMUA                      = NaN(length(params.finalYaxis),nSessions);
% outMeanCSD                      = NaN(length(params.finalYaxis),nTimebins,nSessions);
% outMeanERP                      = NaN(length(params.finalYaxis),nTimebins,nSessions);

%% Get average MUA and CSD across sessions:

for iSes = 1:nSessions
    
    fprintf('Aligning MUA and CSD for session %d/%d\n',iSes,nSessions)
    
    %Get the relevant data for each session individually:
    [templfpData] = MOL_getTempPerSes(sessionData.session_ID(iSes),lfpData);
    
    params.Ylim_all                    	= [min(templfpData.sortedChannelDepth)-25 max(templfpData.sortedChannelDepth)+25];
    params.nChannels                    = length(templfpData.ch);
    
    idx_ch                              = templfpData.sortedChannelDepth>params.chdepthmin-50        & templfpData.sortedChannelDepth<params.chdepthmax+50 & templfpData.sortedisgood;
    idx_out                             = params.finalYaxis>min(templfpData.sortedChannelDepth)-50   & params.finalYaxis<max(templfpData.sortedChannelDepth)+50;
    
    outMeanMUA(idx_out,iSes)            = tointerpol2(templfpData.sortedHF_PWR(idx_ch,:),templfpData.sortedChannelDepth(idx_ch),params.finalYaxis(idx_out));
%     outMeanMUA(idx_out,iSes)            = tointerpol2(templfpData.HF_PWR(idx_ch,:),templfpData.sortedChannelDepth(idx_ch),params.finalYaxis(idx_out));
    
    idx_ch                              = templfpData.sortedChannelDepth>params.chdepthmin-50        & templfpData.sortedChannelDepth<params.chdepthmax+50;
    outMeanCSD(idx_out,:,iSes)          = tointerpol2(templfpData.CheckerCSD(idx_ch,:),templfpData.sortedChannelDepth(idx_ch),params.finalYaxis(idx_out));
    outMeanERP(idx_out,:,iSes)          = tointerpol2(templfpData.CheckerERP(idx_ch,:),templfpData.sortedChannelDepth(idx_ch),params.finalYaxis(idx_out));
    
end

%% Renormalize MUA for each session:

for iSes = 1:nSessions
    outMeanMUA(:,iSes)            = outMeanMUA(:,iSes) - min(outMeanMUA(:,iSes));
    outMeanMUA(:,iSes)            = outMeanMUA(:,iSes) / max(outMeanMUA(:,iSes));
end

%% PLOT AVERAGE MUA AND CSD:

figure; set(gcf,'units','normalized','Position',[0.32 0.45 0.41 0.34],'color','w');

subplot(1,2,1); hold all;

meantoplot_MUA                      = nanmedian(outMeanMUA,2);
% errortoplot_MUA                     = nanstd(outMeanMUA,[],2) / sqrt(nSessions);

for iSes = 1:nSessions
    plot(outMeanMUA(:,iSes),params.finalYaxis,'-','Color',[0.6 0.6 0.6],'LineWidth',0.05)
end

plot(meantoplot_MUA,params.finalYaxis,'k-','LineWidth',2)
% shadedErrorBar(meantoplot_MUA,params.finalYaxis,errortoplot_MUA,{'k-','LineWidth',2},0)
xlabel('Normalized MUA power')
ylim([-1000 0])
set(gca,'YDir','normal','XTick',[0 1])
ylabel('Cortical depth (um)','FontSize', 15)

%Panel 3: CSD power
subplot(1,2,2)

meancsd                             = squeeze(nanmean(outMeanCSD,3)); %get the mean over sessions of this split condition
meanerp                             = squeeze(nanmean(outMeanERP,3));

params.cscale                       = [-max(max(max(meancsd)))*0.92 max(max(max(meancsd)))*0.92];

imagesc(params.xtime,params.finalYaxis,meancsd,params.cscale); hold on;
switch params.colormap_CSD
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

offsetmat = repmat(params.finalYaxis',1,length(params.xtime));
% plot(params.xtime,meanerp*20e4 + offsetmat,'k','LineWidth',0.5); hold on;
idx_subsample = 1:10:params.nTimebins;
plot(params.xtime(idx_subsample),meanerp(:,idx_subsample)*20e4 + offsetmat(:,idx_subsample),'k','LineWidth',0.5); hold on;

plot([0 0],ylim,'k','LineWidth',1);
xlim([params.xtime(1) params.xtime(end)]);
set(gca,'XTick',params.xticks,'XTickLabels',params.xticklabels)
xlabel('Time from stimulus (ms)','FontSize', 15)
ylim([-1000 0])
set(gca,'YDir','normal','XTick',[0 1])

% export_fig(fullfile(params.savedir,'ExSession_%s',tempsessionData.session_ID{1})),'-eps','-nocrop')
export_fig(fullfile(params.savedir,'CSD_MUA_SessionAverage'),'-eps')
% 