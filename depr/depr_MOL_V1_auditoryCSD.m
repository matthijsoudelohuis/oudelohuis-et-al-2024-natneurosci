%% MOL_V1_auditoryCSD
startover

%% Parameters:

global params
%General parameters
params.showFig              = 1;
params.t_pre                = -0.1e6; %All time in microseconds
params.t_post               = .4e6;  %All time in microseconds

params.xticks               = params.t_pre:1e5:params.t_post;
params.xticklabels          = params.xticks*1e-3;

% CSD parameters:
params.flipsignLFP          = 1;    %Whether to flip the LFP because NLX has inversed polarity (so flip again to make comparable to others)
params.csdmethod            = 'NicholsonFreeman'; %'NicholsonFreeman'; %{'NicholsonFreeman' 'DoubleDiff'}
params.colormap             = 'parula'; %redblue or parula
params.sinkhot              = 1;    %Whether sink is shown as red hot colors (2 conventions exist in csd world)
params.conductivity         = 0.4;  %in Siemens m-1

% There are filtering before making the CSD:
% Parameters for Butterworth filter
params.UseButter           = 1;
params.lp_butter           = 100;  %Low pass filter (Hz) (Butter)
params.ord_butter          = 4;   %Butterworth filter order

%Parameters for Kaiser filter
params.UseKaiser           = 0;
params.lp_kaiser           = 80;  %Low pass filter (Hz) (Kaiser)
params.hp_kaiser           = 1;     %High pass filter (Hz) (Kaiser)
params.dp_kaiser           = 1;    %Transition bandwidth (Kaiser only)

params.area                     = {'V1'};

params.nSplits                  = 2;

params.interpolate              = 1;

params.chdepthmin               = -1000;
params.chdepthmax               = 0;
params.resolution               = 20;
params.finalYaxis               = params.chdepthmin:params.resolution:params.chdepthmax;

% params.chdepthmin               = 0;
% params.chdepthmax               = 1000;
% params.resolution               = 20;
% params.finalYaxis               = params.chdepthmin:params.resolution:params.chdepthmax;

params.fs                       = 1024;

params.xtime                    = params.t_pre:1e6/params.fs:params.t_post(1) - 1e6/params.fs;

params.depthcorrection          = 50;

params.Experiments          = {'ChangeDetectionConflict' 'VisOnlyTwolevels' 'ChangeDetectionConflictDecor'};
params.ExperimentLabels     = {'NE' 'UST' 'MST'};

%% Load data:
% [Data]          = MOL_GetData('E:','CHDET',{'ChangeDetectionConflict' 'VisOnlyTwolevels' 'ChangeDetectionConflictDecor'},{},{},{'sessionData' 'trialData' 'lfpData_lazy'});
[Data]          = MOL_GetData('E:','CHDET',params.Experiments,{'2044' '2045'},{},{'sessionData' 'trialData' 'lfpData_lazy'});
sessionData     = Data.sessionData;
lfpData         = Data.lfpData;

%%
filterses                           = unique(lfpData.session_ID(strcmp(lfpData.area,params.area)));

[sessionData,lfpData]               = MOL_getTempPerSes(filterses,sessionData,lfpData);
fprintf('Filtered %d/%d sessions\n',length(sessionData.session_ID),length(Data.sessionData.session_ID));

%% 
params.manualChannelFiltering   = 1;
params.savedir                  = 'E:\Matlab\MOL_Analysis\MOL_AuditoryV1\5LaminarV1CSD\MovArt\';
params.overwrite                = 0;

if params.manualChannelFiltering
    MOL_identifyMovArtefac(sessionData);
end

%%
nSessions                       = length(sessionData.session_ID);

temprunlist = {};

for iSes = 1:nSessions
    load(fullfile(params.savedir,['MovArt_' sessionData.session_ID{iSes} '.mat']),'list')
    temprunlist = [temprunlist; list];
end

lfpfields     = fieldnames(lfpData);
idx             = ismember(lfpData.channel_ID,temprunlist);
for iF = 1:length(lfpfields)
    lfpData.(lfpfields{iF})    = lfpData.(lfpfields{iF})(idx);
end

[sessionData,lfpData]                       = MOL_getTempPerSes(unique(lfpData.session_ID),sessionData,lfpData);
fprintf('Filtered %d/%d sessions\n',length(sessionData.session_ID),nSessions);

%%

params.nSplits                  = 6;

splitlabels                     = {'Vmax' 'Amax' 'Vthr' 'Athr' 'ACmax' 'Amiss'};

%% Init output: interpolated CSD matrix:
nSessions                       = length(sessionData.session_ID);
nTimebins                       = length(params.xtime);

outMeanCSD                      = NaN(length(params.finalYaxis),nTimebins,nSessions,params.nSplits);
outMeanERP                      = NaN(length(params.finalYaxis),nTimebins,nSessions,params.nSplits);

%% For each session compute CSD:

for iSes = 1:nSessions
% for iSes = [13 21 56]
% for iSes = [57 59 60]
    fprintf('Computing CSD for session %d/%d\n',iSes,nSessions)
    
    %Get the relevant data for each session individually with the actual lfp data:
    [Data]              = MOL_GetData('E:','CHDET',{'ChangeDetectionConflict' 'ChangeDetectionConflictDecor' 'VisOnlyTwolevels'},{},sessionData.Rec_datetime(iSes),{'sessionData' 'trialData_newtrials' 'lfpData'});
    tempsessionData     = Data.sessionData;
    temptrialData       = Data.trialData;
    templfpData         = Data.lfpData;
    temptrialData       = MOL_RemoveLastnTrials(temptrialData,5);

    keepidx = strcmp(templfpData.area,params.area);
    lfpfields = fieldnames(templfpData);
    for iF = 1:length(lfpfields)
        templfpData.(lfpfields{iF}) = templfpData.(lfpfields{iF})(keepidx,:);
    end
    
    templfpData.sortedChannelDepth = -templfpData.sortedChannelDepth;
    if all(templfpData.sortedChannelDepth>0)
        templfpData.sortedChannelDepth = -templfpData.sortedChannelDepth;
    end
    
    templfpData.sortedChannelDepth = templfpData.sortedChannelDepth + params.depthcorrection;

    %set channels to bad if manually identified as having movement artefacts:
    normalisbad                     = ~ismember(templfpData.channel_ID,temprunlist);
    idx                             = ismember(templfpData.sortedChannelNum,templfpData.ch(normalisbad));
    templfpData.sortedisgood(idx)   = 0;
%     templfpData.sortedChannelNum(idx)
%     templfpData.lfpgood(templfpData.hasMovArtefacts==1)         = 0;
    
    if ~isempty(templfpData.signal) && any(templfpData.lfpgood)
        
        if params.flipsignLFP
            templfpData.signal = cellfun(@(x) -x, templfpData.signal, 'UniformOutput',false);
        end
        
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
        
        splits          = {};
        splits{1}       = strcmp(temptrialData.trialType,'X') & temptrialData.visualOriChangeNorm==3 & temptrialData.hasphotostim~=1;
        splits{2}       = strcmp(temptrialData.trialType,'Y') & temptrialData.audioOctChangeNorm==3 & temptrialData.hasphotostim~=1;
        
        splits{3}       = strcmp(temptrialData.trialType,'X') & temptrialData.visualOriChangeNorm==2 & temptrialData.hasphotostim~=1;
        splits{4}       = strcmp(temptrialData.trialType,'Y') & temptrialData.audioOctChangeNorm==2 & temptrialData.hasphotostim~=1;
        
        splits{5}       = strcmp(temptrialData.trialType,'Y') & temptrialData.audioOctChangeNorm==3 & temptrialData.vecResponse==3 & temptrialData.hasphotostim~=1;
        splits{6}       = temptrialData.audioOctChangeNorm==3 & temptrialData.hasphotostim~=1;
        
        params.Ylim_all                    	= [min(templfpData.sortedChannelDepth)-25 max(templfpData.sortedChannelDepth)+25];
        params.nChannels                    = length(templfpData.ch);
        
        templfpData.ts                      = repmat(templfpData.t_start(1):1/templfpData.fs(1)*1e6:templfpData.t_end(1),1,1);
        
        templfpData.sortedChannelNum        = mod(templfpData.sortedChannelNum-1,params.nChannels)+1;
        templfpData.sortedChannelidx        = mod(templfpData.sortedChannelidx-1,params.nChannels)+1;

        templfpData.sortedSignal            = cell2mat(templfpData.signal);
        templfpData.sortedSignal            = templfpData.sortedSignal(templfpData.sortedChannelidx,:);
        
        params.ChannelSel                   = 1:params.nChannels;
        
        if strcmp(tempsessionData.session_ID,'20212019070104')
            templfpData.sortedChannelDepth = templfpData.sortedChannelDepth-150;
        elseif strcmp(tempsessionData.session_ID,'20232019070170')
            templfpData.sortedChannelDepth = templfpData.sortedChannelDepth+100;
        elseif strcmp(tempsessionData.session_ID,'20282019121322')
            templfpData.sortedChannelDepth = templfpData.sortedChannelDepth+100;
        elseif strcmp(tempsessionData.session_ID,'20282019121626')
            templfpData.sortedChannelDepth = templfpData.sortedChannelDepth+100;
        elseif strcmp(tempsessionData.session_ID,'20292019121240')
            templfpData.sortedChannelDepth = templfpData.sortedChannelDepth+50;
        elseif strcmp(tempsessionData.session_ID,'20292019121343')
            templfpData.sortedChannelDepth = templfpData.sortedChannelDepth+50;
        end
        
        for iSplit = 1:length(splits)
            events_ts   = temptrialData.stimChange(splits{iSplit});
            
            [meancsd,meanerp] = MOL_CSD(params, events_ts, templfpData);
            
%             idx_ch = templfpData.sortedChannelDepth>params.chdepthmin-25        & templfpData.sortedChannelDepth<params.chdepthmax+25;
%             idx_out = params.finalYaxis>min(templfpData.sortedChannelDepth)-25     & params.finalYaxis<max(templfpData.sortedChannelDepth)+25;
            
            idx_ch = templfpData.sortedChannelDepth>params.chdepthmin-25        & templfpData.sortedChannelDepth<params.chdepthmax+25;
            idx_out = params.finalYaxis>min(templfpData.sortedChannelDepth)-25     & params.finalYaxis<max(templfpData.sortedChannelDepth)+25;
            
            if params.interpolate
                outMeanCSD(idx_out,:,iSes,iSplit)                        = tointerpol2(meancsd(idx_ch,:),templfpData.sortedChannelDepth(idx_ch),params.finalYaxis(idx_out));
                outMeanERP(idx_out,:,iSes,iSplit)                        = tointerpol2(meanerp(idx_ch,:),templfpData.sortedChannelDepth(idx_ch),params.finalYaxis(idx_out));
            end
            
        end
    end
end


%%

filterses  = {     '20202019063070'
    '20202019070175'
    '20202019070381'
    '20202019070483'
    '20212019063001'
    '20212019070104'
    '20212019070413'
    '20222019062636'
    '20222019062739'
    '20222019062841'
    '20232019070170'
    '20232019070376'
    '20232019071684'
    '20232019071787'
    '20262020012103'
    '20302020011448'
    '20302020011551'
    '20302020012261'
    '20312020011491'
    '20312020011594'
    '10082019030758'
    '10082019030869'
    '10082019031271'
    '10082019031373'
    '10092019031200'
    '10112019040515'
    '10122019041022'
    '10122019041126'
    '10132019041036'
    '20322019101795'
    '20332019102302'
    '20332019102405'
    '20352019121876'};

%% Figures and analyses:

%initialize figure
figure; set(gcf,'units','normalized','Position',[0.35   0.3    0.45    0.3],'color','w');

% filterses                           = unique(lfpData.session_ID(lfpData.L4_Sink==1));
% filterses                           = unique(lfpData.session_ID(lfpData.L4_Sink==1 & ismember(lfpData.session_ID,sessionData.session_ID(strcmp(sessionData.Experiment,{'ChangeDetectionConflict'})))));
% filterses                           = unique(lfpData.session_ID(lfpData.L4_Sink==1 & ismember(lfpData.session_ID,sessionData.session_ID(strcmp(sessionData.Experiment,{'ChangeDetectionConflict'})))));
% % filterses                           = sessionData.session_ID(strcmp(sessionData.Experiment,{'ChangeDetectionConflict'}));
% filterses                           = sessionData.session_ID(strcmp(sessionData.Experiment,{'ChangeDetectionConflictDecor'}));
% filterses                           = sessionData.session_ID(strcmp(sessionData.Experiment,{'VisOnlyTwolevels'}));
% filterses                           = unique(lfpData.session_ID(lfpData.L5_maxMUA==1));

idx_ses                             = true(size(sessionData.session_ID));
% idx_ses                             = ismember(sessionData.session_ID,filterses);

%compute mean only to set the colorbarscale:
meancsd                             = squeeze(nanmean(outMeanCSD(:,:,idx_ses,:),3));
params.csdscale                       = [-max(max(max(meancsd)))*0.92 max(max(max(meancsd)))*0.92];
params.csdscale                       = [-max(max(max(meancsd)))*0.9 max(max(max(meancsd)))*0.9];
% params.csdscale                       = [-max(max(max(meancsd)))*0.8 max(max(max(meancsd)))*0.8];
params.erpscale                   = 150e4;

%For each split condition:
for iSplit = 1:2%params.nSplits
%     subplot(1,params.nSplits,iSplit);
    subplot(1,2,iSplit);
    meancsd                             = squeeze(nanmean(outMeanCSD(:,:,idx_ses,iSplit),3)); %get the mean over sessions of this split condition
%     meancsd                             = squeeze(nanmean(abs(outMeanCSD(:,:,idx_ses,iSplit)),3)); %get the mean over sessions of this split condition
    meanerp                             = squeeze(nanmean(outMeanERP(:,:,idx_ses,iSplit),3));
    params.csdscale                       = [-max(max(max(meancsd)))*0.85 max(max(max(meancsd)))*0.85];
    params.csdscale                       = round(params.csdscale,1);
%     params.erpscale                     =  1 / max(max(max(meanerp)))*0.92 * 50;

    %Plot the CSD:
    imagesc(params.xtime,params.finalYaxis,meancsd,params.csdscale); hold on;
    xlabel('Time from stimulus (s)','FontSize', 15)
    set(gca,'YDir','normal')
    colorbar();
    
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
    plot(params.xtime,meanerp*params.erpscale + offsetmat,'k','LineWidth',0.5); hold on; %plot mean erp with offset
    %Figure make up:
    title(splitlabels{iSplit},'FontSize',15)
    ylabel('Cortical Depth','FontSize', 15)
    plot([0 0],ylim,'k','LineWidth',1);
    xlim([params.t_pre params.t_post]);
    set(gca,'XTick',params.xticks,'XTickLabels',params.xticklabels)
    xlabel('Time from stimulus (ms)','FontSize', 15)
end

%% Figure for individual CSD plots:
%initialize figure
% figure; set(gcf,'units','normalized','Position',[-0.81   0.23    1.80    0.6],'color','w');
figure; set(gcf,'units','normalized','Position',[0.35   0.3    0.45    0.55],'color','w');

sesselec                            = 1:10;
% sesselec                            = 11:20;
% sesselec                            = find(ismember(sessionData.session_ID,filterses));
% sesselec                            = [21 39 49 61];
sesselec                            = [39 49];
nSelec                              = length(sesselec);

for iSes = 1:nSelec
    meancsd                             = squeeze(nanmean(outMeanCSD(:,:,sesselec(iSes),:),3)); %get the mean over sessions of this split condition
    params.csdscale                       = [-max(max(max(meancsd)))*0.92 max(max(max(meancsd)))*0.92];
    
    for iSplit = 1:2%params.nSplits%For each split condition:
%         subplot(params.nSplits,nSelec,(iSplit-1)*nSelec+iSes);
        subplot(2,nSelec,(iSplit-1)*nSelec+iSes);
        meancsd                             = squeeze(nanmean(outMeanCSD(:,:,sesselec(iSes),iSplit),3)); %get the mean over sessions of this split condition
        meanerp                             = squeeze(nanmean(outMeanERP(:,:,sesselec(iSes),iSplit),3));
        params.csdscale                     = [-max(max(max(meancsd)))*0.92 max(max(max(meancsd)))*0.92];
        params.csdscale                       = round(params.csdscale,1);
        params.erpscale                     =  1 / max(max(max(meanerp)))*0.92 * 50;
        
        if ~all(isnan(meancsd(:)))
            %Plot the CSD:
            imagesc(params.xtime,params.finalYaxis,meancsd,params.csdscale); hold on;
                set(gca,'YDir','normal')
            colorbar
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
            plot(params.xtime,meanerp*params.erpscale + offsetmat,'k','LineWidth',0.5); hold on; %plot mean erp with offset
            %Figure make up:
            if iSplit==1
%                 title(sessionData.session_ID(iSes),'FontSize',15);
                title(sesselec(iSes),'FontSize',15);
            end
            %         ylabel('Depth','FontSize', 15)
            plot([0 0],ylim,'k','LineWidth',1);
            xlim([params.t_pre params.t_post]);
            set(gca,'XTick',[],'YTick',[])
            set(gca,'XTick',params.xticks,'XTickLabels',params.xticklabels)
        end
    end
end
% tightfig();

%%

params.normalizeCSD = 0;

%initialize figure
csdfig                              = figure; set(gcf,'units','normalized','Position',[0.2 0.2 0.3 0.4],'color','w'); hold all;

idx_time                            = params.xtime>10e3 & params.xtime<70e3;
% idx_time                            = params.xtime>50e3 & params.xtime<100e3;
% idx_time                            = params.xtime>0e3 & params.xtime<200e3;

% idx_ses                             = ismember(sessionData.session_ID,filterses);
idx_ses                             = true(size(sessionData.session_ID));

datatoplot = squeeze(nanmean(abs(outMeanCSD(:,idx_time,:,:)),2)); %get the mean over time

% idx_time_vis                        = params.xtime>40e3 & params.xtime<110e3;
% idx_time_aud                        = params.xtime>10e3 & params.xtime<60e3;
% datatoplot = NaN(length(params.finalYaxis),nSessions,params.nSplits);
% datatoplot(:,:,[1 3]) = squeeze(nanmean(abs(outMeanCSD(:,idx_time_vis,:,[1 3])),2)); %get the mean over time
% datatoplot(:,:,[2 4 5 6]) = squeeze(nanmean(abs(outMeanCSD(:,idx_time_aud,:,[2 4 5 6])),2)); %get the mean over time

if params.normalizeCSD
    for iSes = 1:nSessions
        for iSplit = 1:params.nSplits
            %         datatoplot(:,iSes,:) = datatoplot(:,iSes,:) ./ prctile(tempdata(:),95);
            datatoplot(:,iSes,iSplit) = datatoplot(:,iSes,iSplit) / max(datatoplot(:,iSes,iSplit));
        end
    end
end

params.colors_trialtypes = {[0 0 0.8] [0.8 0.2 0.2]};
params.colors_trialtypes = {[0 0 0.8] [0.8 0.2 0.2] [0.2 0.2 1] [1 0.2 0.2] [1 0 1] [0.4 0.2 0.2]};
% datatoplot = squeeze(nanmean(outMeanCSD(:,idx_time,idx_ses,:),2)); %get the mean over time

% for iSplit = 1:params.nSplits
for iSplit = [1 2]
%     meancsd                             = squeeze(nanmean(outMeanCSD(:,:,idx_ses,iSplit),3)); %get the mean over sessions of this split condition
%     meanerp                             = squeeze(nanmean(outMeanERP(:,:,idx_ses,iSplit),3));

    meantoplot = nanmean(datatoplot(:,idx_ses,iSplit),2);
    
    for iSes = 1:nSessions %plot individual sessions
%         plot(datatoplot(:,iSes,iSplit),params.finalYaxis,'-','LineWidth',0.1,'Color',params.colors_trialtypes{iSplit})
    end
    %Plot the mean time-averaged CSD amplitude:
    plot(meantoplot,params.finalYaxis,'-','LineWidth',2,'Color',params.colors_trialtypes{iSplit})
    xlabel('|CSD|','FontSize', 15)
    set(gca,'YDir','normal') %flip the direction of imagesc
end

params.normalizeCSD = 0;

%% 
csdfig                              = figure; set(gcf,'units','normalized','Position',[0.2 0.2 0.3 0.4],'color','w'); hold all;

idx_time                            = params.xtime>0e3 & params.xtime<50e3;
% idx_time                            = params.xtime>50e3 & params.xtime<100e3;
% idx_time                            = params.xtime>0e3 & params.xtime<200e3;

% idx_ses                             = ismember(sessionData.session_ID,filterses);
idx_ses                             = true(size(sessionData.session_ID));

layerborders = [0 -350 -550 -1000];
layerlabels = {'SG' 'G' 'IG'};

datatoplot = NaN(nSessions,params.nSplits,3);
for iL = 1:3
    idx_depth = params.finalYaxis<=layerborders(iL) & params.finalYaxis>=layerborders(iL+1); 
    datatoplot(:,:,iL) = squeeze(nanmean(nanmean(abs(outMeanCSD(idx_depth,idx_time,:,:)),2),1)); %get the mean over time
    
end

meantoplot = squeeze(nanmean(datatoplot(:,[1 2],:),1));
errortoplot = squeeze(nanstd(datatoplot(:,[1 2],:),[],1));

bar([1 2 3],meantoplot(2,:),'r')
errorbar([1 2 3],meantoplot(2,:),errortoplot(2,:),'k')


%% initialize figure
csdfig                              = figure; set(gcf,'units','normalized','Position',[0.2 0.2 0.25 0.25],'color','w'); hold all;


datatoplot = squeeze(nanmean(abs(outMeanCSD(:,:,:,:)),1)); %get the mean over depth

params.colors_trialtypes = {[0 0 0.8] [0.8 0.2 0.2]};
params.colors_trialtypes = {[0 0 0.8] [0.8 0.2 0.2] [0.2 0.2 1] [1 0.2 0.2] [1 0 1] [0.4 0.2 0.2]};

for iSplit = [1 2]
    meantoplot = nanmean(datatoplot(:,:,iSplit),2);
    %Plot the mean time-averaged CSD amplitude:
    plot(params.xtime(1:end-2),meantoplot(1:end-2),'-','LineWidth',2,'Color',params.colors_trialtypes{iSplit})
end

%Figure make up:
legend(splitlabels(1:2),'FontSize',15); legend boxoff;
ylabel('|CSD|','FontSize', 15)
xlabel('Time from stimulus (ms)','FontSize', 15)
plot([0 0],ylim,'k:','LineWidth',1);
xlim([params.t_pre params.t_post]);
set(gca,'XTick',params.xticks,'XTickLabels',params.xticklabels)
%     xlabel('Time from stimulus (ms)','FontSize', 15)
MOL_prepfigAI


