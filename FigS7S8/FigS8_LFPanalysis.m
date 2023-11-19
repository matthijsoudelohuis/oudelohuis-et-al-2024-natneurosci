%% Oude Lohuis et al. 2023 Nat Neurosci
% Triple dissociation of auditory, visual, and motor processing in primary visual cortex
% MOL (C) 2023

% This script analyzes the LFP data in the audiovisual change detection 
% using ERP analysis and in a laminar fashion (computing CSD)

% The LFP data underlying the analyses is too large to include in the repository
% As well as the data being supplementary to the main analyses

% MOL_V1_auditoryCSD
startover

%% Parameters:

global params
%General parameters
params.showFig              = 1;
params.t_pre                = -0.2e6; %All time in microseconds
params.t_post               = .8e6;  %All time in microseconds

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
params.nTimebins                = length(params.xtime);

params.depthcorrection          = 50;

params.Experiments              = {'ChangeDetectionConflictDecor' 'VisOnlyTwolevels' 'ChangeDetectionConflict'}; %Which versions of the task to load data from
params.ExperimentLabels         = {'NE' 'UST' 'MST'}; %Labels for the different experiments
params.nExperiments             = length(params.Experiments);

params.savedir                  = 'E:\Documents\PhD\Figures\Project CHDET\Results - auV1\5LaminarCircuit\CSD';

%% Load data:
% [Data]          = MOL_GetData('E:','CHDET',{'ChangeDetectionConflict' 'VisOnlyTwolevels' 'ChangeDetectionConflictDecor'},{},{},{'sessionData' 'trialData_newtrials' 'videoData' 'lfpData_lazy'});
% [Data]          = MOL_GetData('E:','CHDET',params.Experiments,{'1008' '1009' '2030' '2031'},{},{'sessionData' 'trialData_newtrials' 'videoData'});
[Data]          = MOL_GetData('E:','CHDET',params.Experiments,{},{},{'sessionData' 'trialData_newtrials' 'videoData'});
sessionData     = Data.sessionData;
trialData       = Data.trialData;
videoData       = Data.videoData;
% lfpData         = Data.lfpData;

%% Remove last n trials:
trialData           = MOL_RemoveLastnTrials(trialData,5);

%% Filter sessions that have lfp in V1:
% filterses                           = unique(lfpData.session_ID(strcmp(lfpData.area,params.area)));
filterses       = unique(sessionData.session_ID(strcmp(sessionData.Probe1_Area,params.area) | strcmp(sessionData.Probe2_Area,params.area) | strcmp(sessionData.Probe3_Area,params.area)));

% [sessionData,trialData,videoData,lfpData]     = MOL_getTempPerSes(filterses,sessionData,trialData,videoData,lfpData);
[sessionData,trialData,videoData]     = MOL_getTempPerSes(filterses,sessionData,trialData,videoData);
fprintf('Filtered %d/%d sessions\n',length(sessionData.session_ID),length(Data.sessionData.session_ID));

% %% Filter sessions that have recordings in V1:
% filterses                           = unique(lfpData.session_ID(strcmp(lfpData.area,params.area)));
% [sessionData,trialData,videoData,lfpData]     = MOL_getTempPerSes(filterses,sessionData,trialData,videoData,lfpData);
% fprintf('Filtered %d/%d sessions\n',length(sessionData.session_ID),length(Data.sessionData.session_ID));

%% Filter out trials with photostimulation in V1:
sesids              = sessionData.session_ID(sessionData.UseOpto & strcmp(sessionData.PhotostimArea,'V1'));
idx                 = trialData.hasphotostim==1 & ismember(trialData.session_ID,sesids);
trialFields         = fieldnames(trialData);
fprintf('Removed %d/%d trials  with photostim in V1\n',sum(idx),length(trialData.session_ID));
for iF = 1:length(trialFields)
    trialData.(trialFields{iF}) = trialData.(trialFields{iF})(~idx);
end

%% Filter out sessions with muscimol:
sesids              = sessionData.session_ID(cellfun(@isempty,sessionData.MuscimolArea));
fprintf('Removed %d/%d sessions  with muscimol manipulations\n',length(sessionData.session_ID)-length(sesids),length(sessionData.session_ID));
[sessionData, trialData, videoData]        = MOL_getTempPerSes(sesids,sessionData,trialData,videoData);

%% Filter sessions based on containing motSVD data:
sesids              = videoData.session_ID(~cellfun(@isempty,videoData.motSVD));
fprintf('Removed %d/%d sessions without motion from video computed\n',length(sessionData.session_ID)-length(sesids),length(sessionData.session_ID));
[sessionData,trialData,videoData] = MOL_getTempPerSes(sesids,sessionData,trialData,videoData);

%% Filter sessions with high-pass filter on LFP:
sesids              = sessionData.session_ID(sessionData.Filter_low<=1);
fprintf('Removed %d/%d sessions with high pass fitler above 1 Hz\n',length(sessionData.session_ID)-length(sesids),length(sessionData.session_ID));
[sessionData,trialData,videoData] = MOL_getTempPerSes(sesids,sessionData,trialData,videoData);

%% Movement artefacts:
params.manualChannelFiltering   = 1;
params.savedir                  = 'E:\Matlab\MOL_Analysis\MOL_AuditoryV1\9LaminarCSD\MovArt\';
params.overwrite                = 0;

if params.manualChannelFiltering
    MOL_identifyMovArtefac(sessionData);
end

%% Load the list of movement artefact channels and remove sessions with only artefacts:
nSessions = length(sessionData.session_ID);

temprunlist         = {};
temprunses          = {};

for iSes = 1:nSessions
    load(fullfile(params.savedir,['MovArt_' sessionData.session_ID{iSes} '.mat']),'list')
    temprunlist     = [temprunlist; list];
    if ~isempty(list)
        temprunses      = [temprunses; sessionData.session_ID(iSes)];
    end
end
[sessionData,trialData,videoData]                       = MOL_getTempPerSes(temprunses,sessionData,trialData,videoData);
fprintf('Remove %d/%d sessions with only movement artefacts\n',length(sessionData.session_ID)-length(temprunses),nSessions);

% lfpfields     = fieldnames(lfpData);
% idx           = ismember(lfpData.channel_ID,temprunlist);
% for iF = 1:length(lfpfields)
%     lfpData.(lfpfields{iF})    = lfpData.(lfpfields{iF})(idx);
% end

% [sessionData,videoData,lfpData]                       = MOL_getTempPerSes(unique(lfpData.session_ID),sessionData,videoData,lfpData);

%% Report dataset:
fprintf('Dataset: %d sessions, %d trials, %d videos \n',length(sessionData.session_ID),length(trialData.session_ID),length(videoData.session_ID));
for iExp = 1:3
    fprintf('%s: %d\n',params.ExperimentLabels{iExp},sum(strcmp(sessionData.Experiment,params.Experiments{iExp})))
end

%%

save('Dataset5.mat','sessionData','trialData','videoData','-v7.3')









%% Get video ME:
params.nSVDs        = 25;
params.AlignOn      = 'stimChange';
params.videofield   = 'motSVD';

[video_hist,video_hist_tot,video_hist_z,params] = MOL_PSTH_video(sessionData,trialData,videoData,params);

%%










%% Give trialData the unique trialID: (correction with _newtrials of older sessions)
trialData.trial_ID      = strcat(trialData.session_ID,cellstr(num2str(trialData.trialNum,'%03.f'))); %Add unique trial ID (sesID + trialnum)

%% 
exampleses = {
    '10092019030898'
    '10132019041242'
    '20282019121626'
    '20352019121876'};


%% Init output: interpolated CSD matrix:
nExampleses = length(exampleses);
params.nSplits                  = 2;
params.labels_splits            = {'Vmax' 'Amax'};

outMeanCSD                      = NaN(length(params.finalYaxis),params.nTimebins,nExampleses,params.nSplits);
outMeanERP                      = NaN(length(params.finalYaxis),params.nTimebins,nExampleses,params.nSplits);

%% For each session compute CSD:

for iSes = 1:nExampleses
    fprintf('Computing CSD for session %d/%d\n',iSes,nExampleses)
    
    sesidx          = strcmp(sessionData.session_ID,exampleses{iSes});
    
    %Get the relevant data for each session individually
    [tempsessionData,temptrialData,tempvideoData]          = MOL_getTempPerSes(sessionData.session_ID(sesidx),sessionData,trialData,videoData);
    %with the actual lfp data:
    [temp]              = MOL_GetData('E:','CHDET',params.Experiments,{},sessionData.Rec_datetime_new(sesidx),{'sessionData' 'lfpData'});
    templfpData         = temp.lfpData;
    
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
        
        splits                              = {};
        splits{1}                           = strcmp(temptrialData.trialType,'X') & temptrialData.visualOriChangeNorm==3;
        splits{2}                           = strcmp(temptrialData.trialType,'Y') & temptrialData.audioFreqChangeNorm==3;

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
            
            [meancsd,meanerp] = MOL_CSD(params, temptrialData.stimChange(splits{iSplit}), templfpData);

            idx_ch          = templfpData.sortedChannelDepth>params.chdepthmin-25        & templfpData.sortedChannelDepth<params.chdepthmax+25;
            idx_out         = params.finalYaxis>min(templfpData.sortedChannelDepth)-25     & params.finalYaxis<max(templfpData.sortedChannelDepth)+25;
                
            if params.interpolate
                outMeanCSD(idx_out,:,iSes,iSplit)                        = tointerpol2(meancsd(idx_ch,:),templfpData.sortedChannelDepth(idx_ch),params.finalYaxis(idx_out));
                outMeanERP(idx_out,:,iSes,iSplit)                        = tointerpol2(meanerp(idx_ch,:),templfpData.sortedChannelDepth(idx_ch),params.finalYaxis(idx_out));
            end
            
        end
    end
end
fprintf( '\n\nDone. \n\n' )

%% Figure for individual CSD plots:
params.colormap                     = 'parula';

for iSes = [4]%:nExampleses
    meancsd                             = squeeze(outMeanCSD(:,:,iSes,2)); %get the mean over sessions of this split condition
    params.csdscale                       = [-max(max(max(meancsd)))*0.95 max(max(max(meancsd)))*0.95];
    params.erpscale                     =  1 / max(max(max(meanerp)))*0.92 * 50;
    params.erpscale                     =  100e4;
    figure; set(gcf,'units','normalized','Position',[0.05   0.3    0.47    0.3],'color','w');

    for iSplit = 1:params.nSplits%For each split condition:
%         subplot(params.nSplits,nSelec,(iSplit-1)*nSelec+iSes);
%         subplot(2,nSelec,(iSplit-1)*nSelec+iSes);
%         subplot(nSelec,params.nSplits,iSplit+(iSes-1)*params.nSplits);
        subplot(1,params.nSplits,iSplit);
        meancsd                             = squeeze(outMeanCSD(:,:,iSes,iSplit)); %get the mean over sessions of this split condition
        meanerp                             = squeeze(outMeanERP(:,:,iSes,iSplit));
%         params.csdscale                     = [-max(max(max(meancsd)))*0.92 max(max(max(meancsd)))*0.92];
%         params.csdscale                     = round(params.csdscale,1);
        
        if ~all(isnan(meancsd(:)) | (meancsd(:)==0)) 
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
                    else
                        colormap(parula);
                    end
            end
            %To plot the ERP over it, compute an offset:
            offsetmat = repmat(params.finalYaxis',1,length(params.xtime));
            plot(params.xtime,meanerp*params.erpscale + offsetmat,'k','LineWidth',0.5); hold on; %plot mean erp with offset
            %Figure make up:
            if iSplit==1
%                 title(sessionData.session_ID(iSes),'FontSize',15);
                title(iSes,'FontSize',15);
            end
            %         ylabel('Depth','FontSize', 15)
            plot([0 0],ylim,'k','LineWidth',1);
            xlim([params.t_pre params.t_post]);
            set(gca,'XTick',[],'YTick',[])
            set(gca,'XTick',params.xticks,'XTickLabels',params.xticklabels)
            xlim([-100e3 400e3])
        end
    end
    
    %Export figure:
    filename = sprintf('CSD_exSes_%s.eps',exampleses{iSes});
%     export_fig(fullfile(params.savedir,filename),gcf);

end
% tightfig();



%% 



%% Set parameters for the different trial types:
params.nSplits                  = 8;
params.labels_splits            = {'Vmax_lw' 'Vmax_hg' 'Amax_lw' 'Amax_hg'};% 'ACmax' 'Amiss'};

%% Init output: interpolated CSD matrix:
nSessions                       = length(sessionData.session_ID);

outMeanCSD                      = NaN(length(params.finalYaxis),params.nTimebins,nSessions,params.nSplits);
outMeanERP                      = NaN(length(params.finalYaxis),params.nTimebins,nSessions,params.nSplits);

%% For each session compute CSD:

params.minNtrials = 3;

for iSes = 1:nSessions
% for iSes = [13 21 56]
% for iSes = [57 59 60]
    fprintf('Computing CSD for session %d/%d\n',iSes,nSessions)
    
    %Get the relevant data for each session individually
    [tempsessionData,temptrialData,tempvideoData]          = MOL_getTempPerSes(sessionData.session_ID(iSes),sessionData,trialData,videoData);
    %with the actual lfp data:
    [temp]              = MOL_GetData('E:','CHDET',params.Experiments,{},sessionData.Rec_datetime_new(iSes),{'sessionData' 'lfpData'});
    templfpData         = temp.lfpData;
    
    
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
        
        temp            = nanmean(video_hist_z(params.xtime_video>0,strcmp(trialData.session_ID,sessionData.session_ID(iSes))),1)';
        
        splits          = {};
        
        tempsplit       = strcmp(temptrialData.trialType,'X') & temptrialData.visualOriChangeNorm==3 & ismember(temptrialData.visualOriPostChangeNorm,[1 2]);
        splits{1}       = tempsplit & temp<0.5 & temp>-0.5;
        splits{2}       = tempsplit & temp>1;
        
        tempsplit       = strcmp(temptrialData.trialType,'X') & temptrialData.visualOriChangeNorm==3 & ismember(temptrialData.visualOriPostChangeNorm,[3 4]);
        splits{3}       = tempsplit & temp<0.5 & temp>-0.5;
        splits{4}       = tempsplit & temp>1;
        
        tempsplit       = strcmp(temptrialData.trialType,'Y') & temptrialData.audioFreqChangeNorm==3 & ismember(temptrialData.audioFreqPostChangeNorm,[1 2]);
        splits{5}       = tempsplit & temp<0.5 & temp>-0.5;
        splits{6}       = tempsplit & temp>1;
        
        tempsplit       = strcmp(temptrialData.trialType,'Y') & temptrialData.audioFreqChangeNorm==3 & ismember(temptrialData.audioFreqPostChangeNorm,[3 4]);
        splits{7}       = tempsplit & temp<0.5 & temp>-0.5;
        splits{8}       = tempsplit & temp>1;
        
        for i = 1:4
%             sum(splits{(i-1)*2+1})
%             sum(splits{(i-1)*2+2})
            if sum(splits{(i-1)*2+1})<params.minNtrials && sum(splits{(i-1)*2+2})<params.minNtrials
                splits{(i-1)*2+1} = false(size(temptrialData.session_ID));
                splits{(i-1)*2+2} = false(size(temptrialData.session_ID));
            end
        end
                
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
            
            [meancsd,meanerp] = MOL_CSD(params, temptrialData.stimChange(splits{iSplit}), templfpData);

            idx_ch          = templfpData.sortedChannelDepth>params.chdepthmin-25        & templfpData.sortedChannelDepth<params.chdepthmax+25;
            idx_out         = params.finalYaxis>min(templfpData.sortedChannelDepth)-25     & params.finalYaxis<max(templfpData.sortedChannelDepth)+25;
                
            if params.interpolate
                outMeanCSD(idx_out,:,iSes,iSplit)                        = tointerpol2(meancsd(idx_ch,:),templfpData.sortedChannelDepth(idx_ch),params.finalYaxis(idx_out));
                outMeanERP(idx_out,:,iSes,iSplit)                        = tointerpol2(meanerp(idx_ch,:),templfpData.sortedChannelDepth(idx_ch),params.finalYaxis(idx_out));
            end
            
        end
    end
end
fprintf( '\n\nDone. \n\n' )


%% temp for backup:
outMeanCSD2 = outMeanCSD;
outMeanERP2 = outMeanERP;

%%
params.nSplits                  = 4;
outMeanCSD                      = NaN(length(params.finalYaxis),params.nTimebins,nSessions,params.nSplits);
outMeanERP                      = NaN(length(params.finalYaxis),params.nTimebins,nSessions,params.nSplits);

outMeanCSD(:,:,:,1) = nanmean(outMeanCSD2(:,:,:,[1 3]),4);
outMeanCSD(:,:,:,2) = nanmean(outMeanCSD2(:,:,:,[2 4]),4);
outMeanCSD(:,:,:,3) = nanmean(outMeanCSD2(:,:,:,[5 7]),4);
outMeanCSD(:,:,:,4) = nanmean(outMeanCSD2(:,:,:,[6 8]),4);

outMeanERP(:,:,:,1) = nanmean(outMeanERP2(:,:,:,[1 3]),4);
outMeanERP(:,:,:,2) = nanmean(outMeanERP2(:,:,:,[2 4]),4);
outMeanERP(:,:,:,3) = nanmean(outMeanERP2(:,:,:,[5 7]),4);
outMeanERP(:,:,:,4) = nanmean(outMeanERP2(:,:,:,[6 8]),4);

%% FIgure on distribution of motion across trials
edges           = -2:0.1:5;
edgeticks       = edges(1:end-1)+0.05;
idx_trials      = ismember(trialData.trialType,{'X' ' Y'});
temp            = nanmean(video_hist_z(params.xtime_video>0,idx_trials),1)';
ncounts         = histcounts(temp,edges, 'Normalization','count');

figure; set(gcf,'units','normalized','Position',[0.13   0.43    0.14    0.24],'color','w'); hold all;
idx = edgeticks<-0.5;
h = bar(edgeticks(idx),ncounts(idx),'k');
h.FaceColor = [0.7 0.7 0.7];
h.EdgeColor = [0.2 0.2 0.2];
idx = edgeticks>-0.5 & edgeticks < 0.5;
h = bar(edgeticks(idx),ncounts(idx),'k');
h.FaceColor = [0.3 0.9 0.3];
h.EdgeColor = [0.2 0.2 0.2];
idx = edgeticks>0.5 & edgeticks < 1;
h = bar(edgeticks(idx),ncounts(idx),'k');
h.FaceColor = [0.7 0.7 0.7];
h.EdgeColor = [0.2 0.2 0.2];
idx = edgeticks>1;
h = bar(edgeticks(idx),ncounts(idx),'k');
h.FaceColor = [0.7 0.2 0.7];
h.EdgeColor = [0.2 0.2 0.2];
ylabel('Counts')
xlabel('Z-scored video ME')

%Export figure:
filename = sprintf('Hist_video_ME_CSDses.eps');
export_fig(fullfile(params.savedir,filename),gcf);

%%

filterses = sessionData.session_ID(~ismember(1:46,[1:9 14 15 20 27:31 40 44]));
% filterses = sessionData.session_ID(~ismember(1:46,[1:7 14 15 17 22 20 27:31 40 44]));

filterses = {
    '20302020011551'
    '20302020011653'
    '20302020012158'
    '20302020012261'
%     '20312020011594'
    '20442021042025'
    '20442021042434'
    '20442021042435'
    '20442021042838'
    '20452021042849'
    '20452021042850'
    '20452021042952'
    '20452021042953'
%     '20452021043055'
    '10082019030758'
    '10082019030869'
%     '10082019031271'
%     '10082019031373'
    '10092019030898'
    '10092019031200'
    '10092019031302'
    '10122019041022'
    '10122019041126'
    '10122019041229'
    '10132019040934'
    '10132019041242'
    '20282019121117'
    '20282019121322'
    '20282019121626'
    '20292019121240'
    '20292019121343'
    '20352019121876'};

filterses =     {'20262020012103'
    '20302020011448'
    '20302020011551'
    '20302020011653'
    '20302020012158'
    '20302020012261'
    '20442021042025'
    '20442021042230'
    '20442021042434'
    '20442021042435'
    '20442021042838'
    '20452021042041'
    '20452021042849'
    '20452021042850'
    '20452021042952'
    '20452021042953'
    '10092019030898'
    '10092019031200'
    '10092019031302'
    '10122019041022'
    '10122019041126'
    '10122019041229'
    '10132019040934'
    '10132019041242'
    '20282019121322'
    '20282019121626'
    '20292019121137'
    '20292019121343'
    '20352019121876'};

%% Figure of CSD and difference:

%initialize figure
figure; set(gcf,'units','normalized','Position',[0.13   0.13    0.85    0.75],'color','w');
idx_ses                             = ismember(sessionData.session_ID,filterses);
% idx_ses                             = ismember(1:46,[28 29]);

%compute mean only to set the colorbarscale:
meancsd                             = squeeze(nanmedian(outMeanCSD(:,:,idx_ses,:),3));
params.csdscale                     = [-max(max(max(meancsd)))*0.9 max(max(max(meancsd)))*0.9];
params.erpscale                     = 150e4;
params.csdscale                     = round(params.csdscale,1);
params.csdscale                     = [-0.6 0.6];

subplotorder = [1 2 4 5 3 6];

params.labels_splits(5:6) = {'Diff' 'Diff'};

for iSplit = 1:params.nSplits+2
    subh = subplot(2,3,subplotorder(iSplit));
    if iSplit==5
        meancsd                             = squeeze(nanmean(outMeanCSD(:,:,idx_ses,2),3)) - squeeze(nanmean(outMeanCSD(:,:,idx_ses,1),3)); %take difference map
        meanerp                             = squeeze(nanmean(outMeanERP(:,:,idx_ses,2),3)) - squeeze(nanmean(outMeanERP(:,:,idx_ses,1),3));
%         params.csdscale                     = [-max(max(max(meancsd)))*0.9 max(max(max(meancsd)))*0.9];
        params.csdscale                     = [-0.4 0.4];
        params.colormap                     = 'redblue';
    elseif iSplit==6
        meancsd                             = squeeze(nanmean(outMeanCSD(:,:,idx_ses,4),3)) - squeeze(nanmean(outMeanCSD(:,:,idx_ses,3),3)); %take difference map
        meanerp                             = squeeze(nanmean(outMeanERP(:,:,idx_ses,4),3)) - squeeze(nanmean(outMeanERP(:,:,idx_ses,3),3));
%         params.csdscale                     = [-max(max(max(meancsd)))*0.9 max(max(max(meancsd)))*0.9];
        params.csdscale                     = [-0.4 0.4];
        params.colormap                     = 'redblue';
    else
        meancsd                             = squeeze(nanmean(outMeanCSD(:,:,idx_ses,iSplit),3)); %get the mean over sessions of this split condition
        meanerp                             = squeeze(nanmean(outMeanERP(:,:,idx_ses,iSplit),3));
        params.colormap                     = 'parula';
    end

    %Plot the CSD:
    imagesc(params.xtime,params.finalYaxis,meancsd,params.csdscale); hold on;
    xlabel('Time from stimulus (s)','FontSize', 15)
    set(gca,'YDir','normal')
    colorbar();
    
    switch params.colormap
        case 'redblue'
            h       = coolwarm + repmat(0.1-abs(linspace(-0.1,.1,256))',1,3);
            h(h>1)  = 1;
            if params.sinkhot
                h = flipud(h);
            end
            colormap(subh,h);
        case 'parula'
            if params.sinkhot
                colormap(subh,flipud(parula));
            else
                colormap(subh,parula);
            end
    end
    %To plot the ERP over it, compute an offset:
    offsetmat = repmat(params.finalYaxis',1,length(params.xtime));
    plot(params.xtime,meanerp*params.erpscale + offsetmat,'k','LineWidth',0.5); hold on; %plot mean erp with offset
    %Figure make up:
    title(params.labels_splits{iSplit},'FontSize',15)
    ylabel('Cortical Depth','FontSize', 15)
    plot([0 0],ylim,'k','LineWidth',1);
    xlim([params.t_pre params.t_post]);
    set(gca,'XTick',params.xticks,'XTickLabels',params.xticklabels)
    xlabel('Time from stimulus (ms)','FontSize', 15)
    xlim([-100e3 600e3])
end

%Export figure:
filename = sprintf('CSD_averagesessions_n46_movementsplit_diff.eps');
export_fig(fullfile(params.savedir,filename),gcf);

%% Figure for individual CSD plots:

sesselec                            = 1:nSessions;
% sesselec                            = 11:20;
% sesselec                            = find(ismember(sessionData.session_ID,filterses));

nSelec                              = length(sesselec);
params.colormap                     = 'parula';

for iSes = 1:nSelec
    meancsd                             = squeeze(nanmean(outMeanCSD(:,:,sesselec(iSes),:),3)); %get the mean over sessions of this split condition
    params.csdscale                       = [-max(max(max(meancsd)))*0.92 max(max(max(meancsd)))*0.92];
%     params.csdscale                       = [-1 1];
    params.erpscale                     =  1 / max(max(max(meanerp)))*0.92 * 50;
    params.erpscale                     =  100e4;
    figure; set(gcf,'units','normalized','Position',[0.05   0.3    0.92    0.4],'color','w');

    for iSplit = 1:params.nSplits%For each split condition:
        subplot(1,params.nSplits,iSplit);
        meancsd                             = squeeze(nanmean(outMeanCSD(:,:,sesselec(iSes),iSplit),3)); %get the mean over sessions of this split condition
        meanerp                             = squeeze(nanmean(outMeanERP(:,:,sesselec(iSes),iSplit),3));
        
        if ~all(isnan(meancsd(:)) | (meancsd(:)==0)) 
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
                    else
                        colormap(parula);
                    end
            end
            %To plot the ERP over it, compute an offset:
            offsetmat = repmat(params.finalYaxis',1,length(params.xtime));
            plot(params.xtime,meanerp*params.erpscale + offsetmat,'k','LineWidth',0.5); hold on; %plot mean erp with offset
            %Figure make up:
            if iSplit==1
                title(sesselec(iSes),'FontSize',15);
            end
            plot([0 0],ylim,'k','LineWidth',1);
            xlim([params.t_pre params.t_post]);
            set(gca,'XTick',[],'YTick',[])
            set(gca,'XTick',params.xticks,'XTickLabels',params.xticklabels)
            xlim([-100e3 400e3])
        end
    end
end


%% GET CI:
data_ci = NaN(4,params.nTimebins,2);
for iSplit = 1:4
    iSplit %#ok<NOPTS>
    for iT = 1:params.nTimebins
        data_ci(iSplit,iT,:) = bootci(1000,@nanmean,datatoplot(iT,idx_ses,iSplit));%,datatoplot(iT,idx_ses,2))
    end
end

%% initialize figure

datatoplot                          = squeeze(nanmean(abs(outMeanERP(:,:,:,:)),1)); %get the mean over depth
%baseline subtraction:
datatoplot                          = datatoplot - repmat(nanmean(datatoplot(params.xtime<0e6,:,:),1),params.nTimebins,1,1);

idx_ses                             = ismember(sessionData.session_ID,filterses);% & idx_exp;

params.colors_trialtypes = {[0.2 0.2 1] [0 0 0.8]  [1 0.2 0.2] [0.8 0 0] [1 0 1] [0.4 0.2 0.2]};

csdfig                              = figure; set(gcf,'units','normalized','Position',[0.2 0.2 0.38 0.43],'color','w'); hold all;

handles = [];
for iSplit = 1:4
    meantoplot      = nanmean(datatoplot(:,idx_ses,iSplit),2);
    errortoplot     = nanstd(datatoplot(:,idx_ses,iSplit),[],2) / sqrt(sum(idx_ses));
    errortoplot = abs(squeeze(data_ci(iSplit,:,:))' - repmat(meantoplot',2,1));
    %Plot the mean time-averaged CSD amplitude:
    h = shadedErrorBar(params.xtime(1:end-2),meantoplot(1:end-2),errortoplot(:,1:end-2),{'-','LineWidth',2,'Color',params.colors_trialtypes{iSplit}},1);
    handles(iSplit) = h.mainLine; delete(h.edge(1:2));
end

data_vis = nanmean(datatoplot(:,idx_ses,[1 2]),3);
data_aud = nanmean(datatoplot(:,idx_ses,[3 4]),3);

for iT = 1:params.nTimebins
    
    p_vis(iT) = signrank(nanmean(data_vis(params.xtime<0,:),1),data_vis(iT,:)) * params.nTimebins;
    p_aud(iT) = signrank(nanmean(data_aud(params.xtime<0,:),1),data_aud(iT,:)) * params.nTimebins;
    
    p_mot_vis(iT) = signrank(datatoplot(iT,idx_ses,1),datatoplot(iT,idx_ses,2)) * params.nTimebins;
    p_mot_aud(iT) = signrank(datatoplot(iT,idx_ses,3),datatoplot(iT,idx_ses,4)) * params.nTimebins;
    
%     xdata = nanmean(datatoplot(params.xtime<0,idx_ses,1)-datatoplot(params.xtime<0,idx_ses,2),1);
%     p_mot_vis(iT) = signrank(xdata,datatoplot(iT,idx_ses,1)-datatoplot(iT,idx_ses,2));
%     xdata = nanmean(datatoplot(params.xtime<0,idx_ses,3)-datatoplot(params.xtime<0,idx_ses,4),1);
%     p_mot_aud(iT) = signrank(xdata,datatoplot(iT,idx_ses,3)-datatoplot(iT,idx_ses,4));
    
%     p_mot_vis(iT) = signrank(datatoplot(iT,idx_ses,1),datatoplot(iT,idx_ses,2)) * 1;
%     p_mot_aud(iT) = signrank(datatoplot(iT,idx_ses,3),datatoplot(iT,idx_ses,4)) * 1;
end

latvis = params.xtime(find(p_vis<0.05,1));
lataud = params.xtime(find(p_aud<0.05,1));
if latvis; plot(latvis,0.5*1e-4,'v','MarkerSize',8,'MarkerFaceColor',params.colors_trialtypes{1},'MarkerEdgeColor',params.colors_trialtypes{1},'LineWidth',2); text(latvis,0.55e-4,sprintf('%2.1f ms',latvis*1e-3)); end
if lataud; plot(lataud,0.5*1e-4,'v','MarkerSize',8,'MarkerFaceColor',params.colors_trialtypes{3},'MarkerEdgeColor',params.colors_trialtypes{3},'LineWidth',2); text(lataud,0.6e-4,sprintf('%2.1f ms',lataud*1e-3));end

latmot_vis = params.xtime(find((data_ci(1,:,1)>data_ci(2,:,2) | data_ci(2,:,1)>data_ci(1,:,2)) & params.xtime>0,1));
latmot_aud = params.xtime(find((data_ci(3,:,1)>data_ci(4,:,2) | data_ci(4,:,1)>data_ci(3,:,2)) & params.xtime>0,1));

% latmot_vis = params.xtime(find(p_mot_vis<0.05,1));
% latmot_aud = params.xtime(find(p_mot_aud<0.05,1));
if latmot_vis; plot(latmot_vis,0.7*1e-4,'v','MarkerSize',8,'MarkerFaceColor','w','MarkerEdgeColor',params.colors_trialtypes{2},'LineWidth',2); text(latmot_vis,0.65e-4,sprintf('%2.1f ms',latmot_vis*1e-3)); end
if latmot_aud; plot(latmot_aud,0.7*1e-4,'v','MarkerSize',8,'MarkerFaceColor','w','MarkerEdgeColor',params.colors_trialtypes{4},'LineWidth',2); text(latmot_aud,0.7e-4,sprintf('%2.1f ms',latmot_aud*1e-3)); end

%Figure make up:
ylabel('Voltage (\muV)','FontSize', 15)
xlabel('Time from stimulus change (ms)','FontSize', 15)
plot([0 0],ylim,'k:','LineWidth',1);
xlim([params.t_pre params.t_post]);
ylim([-12e-6 70e-6])
set(gca,'XTick',params.xticks,'XTickLabels',params.xticklabels,'YTick',(-1:0.01:6)*1e-3,'YTickLabels',(-1:0.01:6)*1e3)
legend(handles,params.labels_splits,'FontSize',15); legend boxoff;
grid on;
MOL_prepfigAI

