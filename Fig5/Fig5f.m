%% Oude Lohuis et al. 2023 Nat Neurosci
% Triple dissociation of auditory, visual, and motor processing in primary visual cortex
% MOL (C) 2023

% Script analyzes firing rate in V1 with and without muscimol in AC

startover

%% Parameter settings:
params.Experiments          = 'ChangeDetectionConflict'; %Which versions of the task to load data from

params                      = MOL_getColors_CHDET(params);
params.areas                = 'V1';

% histogram parameters:
params                      = params_histresponse_auV1(params);
params.AlignOn              = 'stimChange';         %On which timestamp to align trials to

params.binsize              = 10e3;
params.t_pre                = -1e6;

% Construct bin edges and time axis
params.edges                = [params.t_pre:params.binsize:params.t_post] - params.binsize/2;                    %#ok<NBRAK> Define edges of the bins
params.xtime                = params.t_pre:params.binsize:params.t_post-params.binsize;                        %Define time axis
params.nTimebins            = length(params.xtime); %number of time bins

params.savedir              = 'E:\Documents\PhD\Figures\Project CHDET\Results - auV1\7MuscimolA1\V1responses';

%% Get data:
[Data] = MOL_GetData('E:','CHDET',params.Experiments,{},[],{'sessionData' 'trialData_newtrials' 'spikeData'});
sessionData         = Data.sessionData;
trialData           = Data.trialData;
spikeData           = Data.spikeData;

%% Remove last 20 trials:
trialData           = MOL_RemoveLastnTrials(trialData,20);

%% Remove sessions that are passive: 
sesids                      = sessionData.session_ID(strcmp(sessionData.State,'Behaving'));
fprintf('Removed %d/%d sessions with passive stimulation\n',length(sessionData.session_ID)-length(sesids),length(sessionData.session_ID));
[sessionData, trialData, spikeData]        = MOL_getTempPerSes(sesids,sessionData,trialData,spikeData);

%% Filter out trials with photostimulation in V1:
sesids              = sessionData.session_ID(sessionData.UseOpto & strcmp(sessionData.PhotostimArea,'V1'));
idx                 = trialData.hasphotostim==1 & ismember(trialData.session_ID,sesids);
trialFields         = fieldnames(trialData);
fprintf('Removed %d/%d trials  with photostim in V1\n',sum(idx),length(trialData.session_ID));
for iF = 1:length(trialFields)
    trialData.(trialFields{iF}) = trialData.(trialFields{iF})(~idx);
end

%% Filter neurons on area:
spikeFields     = fieldnames(spikeData);
idx             = ismember(spikeData.area,params.areas);
fprintf('Filtered %d/%d neurons based on area\n',sum(idx),length(spikeData.session_ID));
for iF = 1:length(spikeFields)
    spikeData.(spikeFields{iF}) = spikeData.(spikeFields{iF})(idx);
end

%% Get only sessions with neurons recorded 
[sessionData,trialData,spikeData] = MOL_getTempPerSes(unique(spikeData.session_ID),sessionData,trialData,spikeData);

%% Save dataset:
save('Dataset5_3.mat','params','sessionData','trialData','spikeData')

%% Or start script from saved dataset:
load Dataset5_3.mat

%% Report dataset:
fprintf('Dataset: %d sessions, %d trials, %d neurons\n',length(sessionData.session_ID),length(trialData.session_ID),length(spikeData.session_ID));


%% For each neuron compute firing rate:
params.colors_trialtypes    = [params.colors_visual_opto(1:2) params.colors_audio_opto(1:2)];
params.labels_trialtypes    = {'Vthr' 'Vmax' 'Athr' 'Amax'};

params.nSplits              = 4;
params.minTrialCond         = 5;

params.twin_baseline_start = -0.25e6;
params.twin_baseline_stop   = 0;

nNeurons                    = length(spikeData.ts);
lastsesid                   = []; %memory var to keep track of whether neuron comes from different session and design matrix needs to be reconstructed
fprintf('Computing average Z-scored response for neuron        \n');
zmat                        = NaN(nNeurons,params.nTimebins,params.nSplits);

for iNeuron = 1:nNeurons %Loop over all neurons:
    fprintf(repmat('\b', 1, numel([num2str(iNeuron-1) num2str(nNeurons)])+2));
    fprintf('%d/%d\n',iNeuron,nNeurons);
    
    if ~strcmp(lastsesid,spikeData.session_ID(iNeuron)) %construct new predictor matrix if neuron comes from a new session:
        %Get the relevant data for each session individually:
        temptrialData        = MOL_getTempPerSes(spikeData.session_ID(iNeuron),trialData);
        lastsesid            = spikeData.session_ID(iNeuron); %save this session_ID
    end
    
    %Compute histogram:
    events_ts               = temptrialData.(params.AlignOn);
    hist_mat                = calc_psth(events_ts,spikeData.ts{iNeuron},params);    %Construct histogram matrix:
    
    splits                  = {};
    splits{1}               = ismember(temptrialData.trialType,{'X'}) & temptrialData.visualOriChangeNorm==2;
    splits{2}               = ismember(temptrialData.trialType,{'X'}) & temptrialData.visualOriChangeNorm==3;
    splits{3}               = ismember(temptrialData.trialType,{'Y'}) & temptrialData.audioFreqChangeNorm==2;
    splits{4}               = ismember(temptrialData.trialType,{'Y'}) & temptrialData.audioFreqChangeNorm==3;
    
    for iSplit = 1:params.nSplits %Store the mean response for each of these splits
        if sum(splits{iSplit})>=params.minTrialCond
            zmat(iNeuron,:,iSplit) = nanmean(hist_mat(splits{iSplit},:),1);
        end
    end
end

%% Figure settings:

params.labels_ztrials       = {'Vthr' 'Vmax' 'Athr' 'Amax'};
params.lines_ztrials        = {'-' '-' '-' '-'};

params.lines_mans           = {'-' ':'};

params.colors_ztrials       = [params.colors_ztrials(1:2) {[0.9 0.3 0.1] [0.5 0 0.1]}];

%% Construct indices for neurons:
idx_mus                     = ismember(spikeData.session_ID,sessionData.session_ID(strcmp(sessionData.MuscimolArea,'A1')));

idx_ctr                     = ismember(spikeData.session_ID,sessionData.session_ID(~strcmp(sessionData.MuscimolArea,'A1')));

idx_all                     = [idx_ctr idx_mus];
labels                      = {'Ctrl' 'Mus'};

%% Bootstrap mean and error + stat test of difference:

nSubsampling        = sum(idx_mus)*2;

nBoots              = 10000;

bootmean            = NaN(2,4,params.nTimebins);
booterror           = NaN(2,4,params.nTimebins,2);
bootsamp            = NaN(2,nBoots,params.nTimebins);
bootstat            = NaN(4,params.nTimebins);

bootbaseline        = NaN(4,1);

for iSplit = 1:4 %for each trialtype
    for iMan = 1:2 %for each manipulation / treatment
        
        for iBoot = 1:nBoots
            temp                            = find(idx_all(:,iMan));
            bootsamp(iMan,iBoot,:)          = nanmean(zmat(temp(randi(length(temp),nSubsampling,1)),:,iSplit),1);
        end
        
        bootmean(iMan,iSplit,:)       = prctile(bootsamp(iMan,:,:),50,2);
        booterror(iMan,iSplit,:,:)    = squeeze(prctile(bootsamp(iMan,:,:),[97.5 2.5],2))';
    end
    
    bootbaseline(iSplit)       = prctile(reshape(bootsamp(:,:,params.xtime<0),2*nBoots*sum(params.xtime<0),1),95);
    
    tempdiffmat = squeeze(diff(bootsamp,[],1));
    bootstat(iSplit,:) = sum(tempdiffmat>0,1)/nBoots;
    bootstat(iSplit,bootstat(iSplit,:)>0.5) = 1 - bootstat(iSplit,bootstat(iSplit,:)>0.5);
end

bootstat(bootstat==0) = 1/nBoots; %set to smallest possible significance value

%% Make figure of the mean:

figure; set(gcf,'units','normalized','Position',[0.4 0.41 0.23 0.19],'color','w')
tbl = table();
tbl.Time = params.xtime';

for iSplit = [2 4]
    subplot(1,2,iSplit/2); hold all;
    handles = [];
    for iMan = 1:2
        
        meantoplot = nanmean(zmat(idx_all(:,iMan),:,iSplit),1);
        errortoplot = nanstd(zmat(idx_all(:,iMan),:,iSplit),1)/sqrt(size(zmat(idx_all(:,iMan),:,iSplit),1));
        plot([0 0], [-2 5],'k:','LineWidth',1);
        
        tbl.(sprintf('%s_%s_mean',params.labels_ztrials{iSplit},labels{iMan})) = meantoplot';
        tbl.(sprintf('%s_%s_error',params.labels_ztrials{iSplit},labels{iMan})) = errortoplot';

        if ~all(isnan(meantoplot))
            h = shadedErrorBar(params.xtime,meantoplot,errortoplot,{'k','markerfacecolor',params.colors_ztrials{iSplit-iMan+1},'LineWidth',1,'LineStyle',params.lines_mans{iMan}},0);
            h.mainLine.Color = params.colors_ztrials{iSplit-iMan+1};    h.patch.FaceColor = params.colors_ztrials{iSplit-iMan+1}; delete(h.edge(:));
            handles(end+1) = h.mainLine; hold all; %#ok<SAGROW>
        end
    end
    
    %Statistics:
    markerColors = -log10(bootstat(iSplit,:));
    markerColors(1) = 0;
    scatter(params.xtime,ones(params.nTimebins,1)*-0.1,10,markerColors,'s','fill')
    
    set(gca, 'XTick', -3e6:0.1e6:3e6, 'XTickLabels',(-3e6:0.1e6:3e6)/1e3,'FontSize', 10)
    set(gca, 'YTick', [0 0.5 1 1.5 2], 'FontSize', 10)
    xlim([-0.15e6 0.5e6]);
    ylim([-0.2 2])
    
    title(labels{iMan})
    xlabel('Time (ms)','FontSize', 12)

end

MOL_prepfigAI();
pause(0.1)

export_fig(fullfile(params.savedir,'Zscore_rate_bt_diff'),'-eps','-nocrop')

tempf = figure; 
imagesc(randn(50));
caxis([0 1])
c = colorbar();
c.Ticks = -log10([0.1 0.01 0.001 1/nBoots]) / -log10([1/nBoots]);
c.TickLabels = [0.1 0.01 0.001 1/nBoots];

export_fig(fullfile(params.savedir,'colorbar_musc_sig'),'-eps','-nocrop')
close(tempf)

writetable(tbl,'SourceData_Fig5f_V1_spiking_Muscimol.xlsx')


%% Make figure of the close up of the initial part:

%parameters for zoom in:
xmin = -0.04e6;
xmax = 0.18e6;
xresol = 40e3;
ymin = -0.1;
ymax = 0.3;
yresol = 0.1;

figure; set(gcf,'units','normalized','Position',[0.4 0.41 0.23 0.19],'color','w')

ythr = nanmean(bootbaseline);

crossvals = [];

for iSplit = 1:4
    for iMan = 1:2
        crossvals(iMan,iSplit,1) = params.xtime(find(squeeze(bootmean(iMan,iSplit,:)>ythr)' & params.xtime>0,1));
        crossvals(iMan,iSplit,2) = params.xtime(find(squeeze(booterror(iMan,iSplit,:,1)>ythr)' & params.xtime>0,1));
        crossvals(iMan,iSplit,3) = params.xtime(find(squeeze(booterror(iMan,iSplit,:,2)>ythr)' & params.xtime>0,1));
        
    end
end

fprintf('Onset latency was not significantly different for visual stimuli (Ctrl: %2.0f ms (%2.0f - %2.0f) ms versus Mus: %2.0f ms (%2.0f - %2.0f) ms; Median and bootstrapped 95%%CI; maximal saliency)\n ',crossvals(1,2,:)*1e-3,crossvals(2,2,:)*1e-3)
fprintf('but delayed for auditory stimuli (Ctrl: %2.0f ms (%2.0f - %2.0f) ms versus Mus: %2.0f ms (%2.0f - %2.0f) ms)\n ',crossvals(1,4,:)*1e-3,crossvals(2,4,:)*1e-3)

crossvals(:,:,2) =  crossvals(:,:,2) - crossvals(:,:,1);
crossvals(:,:,3) =  crossvals(:,:,3) - crossvals(:,:,1);

for iSplit = [2 4]
    subplot(1,2,iSplit/2); hold all;
    handles = [];
    
    for iMan = 1:2
        plot([0 0], [-2 5],'k:','LineWidth',1);
        
        if ~all(isnan(bootmean))
            ploterror = abs(squeeze(booterror(iMan,iSplit,:,:) - repmat(bootmean(iMan,iSplit,:),1,1,1,2)));    %from percentiles make error relative to the median:
            plotmedian = squeeze(bootmean(iMan,iSplit,:));    %median
            h = shadedErrorBar(params.xtime,plotmedian,ploterror,{'k','markerfacecolor',params.colors_ztrials{iSplit-iMan+1},'LineWidth',1,'LineStyle',params.lines_mans{iMan}},0);
            
            h.mainLine.Color = params.colors_ztrials{iSplit-iMan+1};    h.patch.FaceColor = params.colors_ztrials{iSplit-iMan+1}; delete(h.edge(:));
            handles(end+1) = h.mainLine; hold all; %#ok<SAGROW>
        end
        errorbar(crossvals(iMan,iSplit,1),ythr+0.01*iMan,crossvals(iMan,iSplit,2),crossvals(iMan,iSplit,3),'Horizontal','Color',params.colors_ztrials{iSplit-iMan+1},'LineWidth',2)
        plot(crossvals(iMan,iSplit,1),ythr+0.01*iMan,'k.','MarkerSize',15)

    end
    
    plot([xmin xmax],[ythr ythr],':','Color',[0.4 0.4 0.4],'LineWidth',0.25)
    xlim([xmin xmax])
    ylim([ymin ymax])
    set(gca,'XTick',xmin:xresol:xmax,'XTickLabels',(xmin:xresol:xmax)*1e-3,'YTick',ymin:yresol:ymax,'YTickLabels',(ymin:yresol:ymax),'XTickLabelRotation',45)

    %Figure make up:
    xlabel('Time (ms)')
    ylabel('Firing rate (z-score)')

    if iMan==1
        title(labels{iMan})
    end
    xlabel('Time (ms)','FontSize', 12)
end

filename = sprintf('PopulationLatency_V1_AV_Musc.eps');
% filename = sprintf('PopulationLatency_V1_AV_Musc.svg');
export_fig(fullfile(params.savedir,filename),gcf);
% saveas(gcf,fullfile(params.savedir,filename));


