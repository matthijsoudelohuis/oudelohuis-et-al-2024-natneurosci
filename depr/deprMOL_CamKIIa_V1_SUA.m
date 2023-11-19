%% 
startover
    
%% Parameters:

params                          = params_histresponse;
params.zscore                   = 0;

%General parameters
params.t_pre                    = -.1e6; %All time in microseconds
params.t_post                   = .1e6;  %All time in microseconds

params.conv_sigma               = 5e3;
params.twin_baseline_start      = -0.2e6;
params.twin_baseline_stop       = 0;

params.xticks                   = params.t_pre:20e3:params.t_post;
params.xticklabels              = params.xticks*1e-3;

% params.yticks                   = [-1000:200:-0];
params.yticks                   = 0:200:1000;
params.yticklabels              = params.yticks;

%Construct bin edges and time axis
params.edges                    = [params.t_pre:params.binsize:params.t_post] - params.binsize/2;                    %#ok<NBRAK> Define edges of the bins
params.xtime                    = params.t_pre:params.binsize:params.t_post-params.binsize;                        %Define time axis
params.nTimebins                = length(params.xtime); %number of time bins

params.area                     = 'V1';

% params.interpolate              = 1;
params.chdepthmin               = 0;
params.chdepthmax               = 1000;

params.AlignOn                  = 'photostimStart';

% params.chdepthmin               = -1000;
% params.chdepthmax               = 200;
params.experiment               = 'OptoOnly';

params.resolution               = 20;
params.finalYaxis               = params.chdepthmin:params.resolution:params.chdepthmax;

% params.fs                       = 1024;

% params.xtime                    = params.t_pre:1e6/params.fs:params.t_post(1) - 1e6/params.fs;
% nTimebins                       = length(params.xtime);

params.exsession = '2019-03-06_12-04-21';
params.exsession = '2019-03-06_12-04-22';
params.exsession = '2019-03-06_12-04-24';
% params.exsession = '2019-03-06_12-04-21';
% params.exsession = '2019-03-07_12-46-60';
% params.exsession = '2019-03-07_10-23-45';
% params.exsession  = '2019-03-06_14-12-20';
% params.exsession  = '2019-04-02_14-34-15';
% params.exsession = '2019-03-07_12-46-57';


%%
params.experiment = 'CheckerboardReversals';
params.exsession = '2019-03-06_12-04-15';
% params.exsession = '2019-03-07_10-23-36';
% params.exsession = '2019-03-07_12-46-53';
params.AlignOn = 'stimStart';

%% 
params.animals                  = {'1008' '1009' '1011'};
params.experiment               = 'OptoOnly';
params.AlignOn                  = 'photostimStart';
% params.animals = {'1008'};

%% Load data:
[Data]          = MOL_GetData('E:','CHDET',params.experiment,params.animals,[],{'sessionData' 'trialData' 'spikeData'});
sessionData     = Data.sessionData;
trialData       = Data.trialData;
spikeData       = Data.spikeData;

%% Load data:
% [Data]          = MOL_GetData('E:','CHDET',params.experiment,[],params.exsession,{'sessionData' 'trialData' 'spikeData'});
% sessionData     = Data.sessionData;
% trialData       = Data.trialData;
% spikeData       = Data.spikeData;


%% Filter neurons on area:
spikeFields     = fieldnames(spikeData);
idx             = ismember(spikeData.area,params.area);
for iF = 1:length(spikeFields)
    spikeData.(spikeFields{iF}) = spikeData.(spikeFields{iF})(idx);
end
fprintf('Filtered %d neurons based on area\n',sum(idx));

%% Compute average response to opto stimulation to see laminar position of response:
nNeurons                = length(spikeData.ts);
lastsesid               = []; %memory var to keep track of whether neuron comes from different session and design matrix needs to be reconstructed
fprintf('Computing average Z-scored response for neuron        \n');
meanmat_opto             = NaN(nNeurons,params.nTimebins);

idx_bsl = params.xtime>-0.02e6 & params.xtime<0;

for iNeuron = 1:nNeurons %Loop over all neurons:
    fprintf(repmat('\b', 1, numel([num2str(iNeuron-1) num2str(nNeurons)])+2));
    fprintf('%d/%d\n',iNeuron,nNeurons);
    
    if ~strcmp(lastsesid,spikeData.session_ID(iNeuron)) %construct new predictor matrix if neuron comes from a new session:
        %Get the relevant data for each session individually:
        temptrialData        = MOL_getTempPerSes(spikeData.session_ID(iNeuron),trialData);
        lastsesid            = spikeData.session_ID(iNeuron); %save this session_ID
    end
    
    %Compute histogram:
    events_ts                   = temptrialData.(params.AlignOn);
    hist_mat                    = calc_psth(events_ts,spikeData.ts{iNeuron},params);    %Construct histogram matrix:
    hist_mat                    = hist_mat - repmat(nanmean(hist_mat(:,idx_bsl),2),1,params.nTimebins);
    meanmat_opto(iNeuron,:)     = mean(hist_mat(:,:),1);
end

%% Make SUA heatmap over laminar depth: (auditory trials)
%Histogram with binning on depth:
binedges                = 0:25:1150;
nBins                   = length(binedges)-1;
bincenters              = binedges(1:end-1)-25/2;

laminardepthfig = figure; set(gcf,'units','normalized','Position',[0.5 0.45 0.2 0.3],'color','w');
hold all; set(laminardepthfig, 'DefaultLineLineWidth', 2);
title('Opto SUA','FontSize',20);


idx_ses                     = sessionData.Photostimpower>=1 & strcmp(sessionData.PhotostimArea,'A1');
idx_spike_ses               = ismember(spikeData.session_ID,sessionData.session_ID(idx_ses));

SUAdepthmap                 = zeros(nBins,params.nTimebins); %Init variable

for iBin = 1:nBins
    idx = idx_spike_ses & spikeData.ChannelY>=binedges(iBin) & spikeData.ChannelY<binedges(iBin+1);
%     idx = spikeData.ChannelY>=binedges(iBin) & spikeData.ChannelY<binedges(iBin+1);
    sum(idx)
    SUAdepthmap(iBin,:) = nanmean(meanmat_opto(idx,:),1);
end

SUAdepthmap(isnan(SUAdepthmap)) = 0;

% Filter temporally and spatially:
spat_filter     = fspecial('gaussian',[50 50],1.1);
spat_filter     = fspecial('gaussian',[25 25],1.1);
SUAdepthmap     = conv2(SUAdepthmap,spat_filter,'same');

params.cscale   = [-max(max(max(SUAdepthmap)))*0.92 max(max(max(SUAdepthmap)))*0.92];

% imagesc(flipud(SUAdepthmap));
imagesc(params.xtime,bincenters,flipud(SUAdepthmap),params.cscale);
set(gca,'YDir','reverse') %flip the direction of imagesc

%Figure make up:
ylabel('Depth','FontSize', 15)
plot([0 0],ylim,'k','LineWidth',1);
xlim([params.t_pre params.t_post]);
xlim([-50e3 params.t_post]);
set(gca,'XTick',params.xticks,'XTickLabels',params.xticklabels,'FontSize',15)
set(gca,'YTick',params.yticks,'YTickLabels',params.yticklabels,'FontSize',15)
xlabel('Time from laser onset (ms)','FontSize', 15)
ylim([params.chdepthmin params.chdepthmax])

%%
ylabel('Depth from dura (um)','FontSize', 15)
xlabel('Time (ms)','FontSize',15)
timeticks = [-0.1e6 0 0.1e6 0.2e6];
idx = find(ismember(params.xtime,timeticks));
% set(gca,'XTick',idx,'XTickLabel',timeticks*1e-3,'fontsize',15,'FontWeight','bold')
set(gca,'XTick',params.xticks,'XTickLabel',params.xticklabels,'fontsize',15,'FontWeight','bold')

% binedges = binedges(1:end-1) + (binedges(2)-binedges(1))/2;
idx = find(ismember(fliplr(binedges(1:end-1)),binticks));
set(gca,'YTick',idx,'YTickLabel',binticks,'fontsize',15,'FontWeight','bold')
set(gca,'fontsize',15,'FontWeight','bold')
