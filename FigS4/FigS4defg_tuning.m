%% Oude Lohuis et al. 2023
startover

%% Parameter settings:
params                      = params_histresponse_auV1; %Parameters for PSTH (All time is in microseconds)
params                      = MOL_getColors_CHDET(params);

params.Experiments          = {'DetectionConflict' }; %Which versions of the task to load data from

params.AlignOn              = 'stimStart';      %On which timestamp to align as t=0

params.trialcategories      = 'OriFreq';

params.allOris              = 0:45:315;
params.allFreqs             = 8e3:1e3:15e3;

params.minNtrials           = 5;

params.alpha                = 0.05;

params.colors_ztrials       = jet(16);

params.labels_ztrials       = num2cell([0:45:315 8:15]);

params.savedir              = 'E:\Documents\PhD\Figures\Project CHDET\Results - auV1\11DetectionTask\Tuning\';

%% Load the data:
[Data] = MOL_GetData('E:','DET',{'DetectionConflict'},{},[],{'sessionData' 'trialData' 'spikeData'});
sessionData     = Data.sessionData;
trialData       = Data.trialData;
spikeData       = Data.spikeData;

%% Get only sessions with neurons recorded:
[sessionData,trialData,spikeData] = MOL_getTempPerSes(unique(spikeData.session_ID),sessionData,trialData,spikeData);
fprintf('Dataset: %d sessions, %d trials, %d neurons\n',length(sessionData.session_ID),length(trialData.session_ID),length(spikeData.session_ID));

%% Save dataset:
save('DatasetS4_3.mat','params','sessionData','trialData','spikeData')

%% Or start script from saved dataset:
load DatasetS4_3.mat
rng(500)

%% Fig S4d,e  - show raster plot of example neurons 
cell_IDs            = {};

%Orientation selective example neurons:
cell_IDs{end+1}     = '3032017100321028'; %very nice

%Frequency selective neurons:
cell_IDs{end+1}     = '3032017100521132'; %neurons from same session, similarly tuned
cell_IDs{end+1}     = '3032017100521167';

params.exportfig            = 0;

params.zscore               = 0;

MOL_plotTuning_OriFreqTuning(sessionData,trialData,spikeData,cell_IDs,params)

%% Main loop to get psth matrix:
params.zscore                   = 1;

nNeurons                        = length(spikeData.ts);

nOris                           = length(params.allOris);
nFreqs                          = length(params.allFreqs);

respmat_ori                     = NaN(nNeurons,nOris);
respmat_ori_trials              = NaN(nNeurons,nOris,100);

respmat_freq                 	= NaN(nNeurons,nFreqs);
respmat_freq_trials             = NaN(nNeurons,nFreqs,100);

lastsesid               = []; %memory var to keep track of whether neuron comes from different session and design matrix needs to be reconstructed
fprintf('Computing average firing rate response for neuron        \n');

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
    
    resp_mat                = calc_resp_from_psth(hist_mat,params);
    
    for iOri = 1:nOris %Store the mean response for each of these splits
        idx                                             = strcmp(temptrialData.trialType,'V') & temptrialData.visualOri==params.allOris(iOri);
        respmat_ori(iNeuron,iOri)                       = mean(resp_mat(idx));
        respmat_ori_trials(iNeuron,iOri,1:sum(idx))     = resp_mat(idx);
    end
    
    for iFreq = 1:nFreqs %Store the mean response for each of these splits
        idx                                             = strcmp(temptrialData.trialType,'A') & temptrialData.audioFreq==params.allFreqs(iFreq);
        respmat_freq(iNeuron,iFreq)                     = mean(resp_mat(idx));
        respmat_freq_trials(iNeuron,iFreq,1:sum(idx))   = resp_mat(idx);
    end
end

%% Get significantly tuned neurons:
%permutation test on shuffled responses:
nShufs              = 1000;
gOSI_ori            = nan(nNeurons,1);
gOSI_ori_shuf       = nan(nNeurons,nShufs);

gOSI_freq           = nan(nNeurons,1);
gOSI_freq_shuf      = nan(nNeurons,nShufs);

fprintf('Computing tuning for neuron        \n');

for iNeuron = 1:nNeurons %Loop over all neurons:
    fprintf(repmat('\b', 1, numel([num2str(iNeuron-1) num2str(nNeurons)])+2));
    fprintf('%d/%d\n',iNeuron,nNeurons);
    
    gOSI_ori(iNeuron,1)     = 1 - circ_var((pi/180) * params.allOris, respmat_ori(iNeuron,:),[],2);
    gOSI_freq(iNeuron,1)    = 1 - circ_var((pi/180) * params.allOris, respmat_freq(iNeuron,:),[],2);
    
    tempresp_ori            = squeeze(respmat_ori_trials(iNeuron,:,:));
    tempresp_freq           = squeeze(respmat_freq_trials(iNeuron,:,:));
    
    for iShuf = 1:nShufs %loop over shuffles
        X               = randperm(numel(tempresp_ori));
        ShuffledData    = reshape(tempresp_ori(X),size(tempresp_ori));
        gOSI_ori_shuf(iNeuron,iShuf) = 1 - circ_var((pi/180) * params.allOris, nanmean(ShuffledData,2)',[],2);
        
        X               = randperm(numel(tempresp_freq));
        ShuffledData    = reshape(tempresp_freq(X),size(tempresp_freq));
        gOSI_freq_shuf(iNeuron,iShuf) = 1 - circ_var((pi/180) * params.allOris, nanmean(ShuffledData,2)',[],2);
    end
end

%% Construct binary index of significantly tuned neurons:
gOSI_ori_sign = gOSI_ori > prctile(gOSI_ori_shuf,(1-params.alpha)*100,2);
gOSI_freq_sign = gOSI_freq > prctile(gOSI_freq_shuf,(1-params.alpha)*100,2);

gOSI_ori_p = sum(gOSI_ori>gOSI_ori_shuf,2)/nShufs;
gOSI_freq_p = sum(gOSI_freq>gOSI_freq_shuf,2)/nShufs;

%% Fig S4f - Show for an example session how the neurons in that session are tuned:
params.exampleSes = {'3032017100541'};
nSessions = length(sessionData.session_ID);

for iSes = find(ismember(sessionData.session_ID,params.exampleSes))'
    idx             = strcmp(spikeData.session_ID,sessionData.session_ID(iSes));
    nSesNeurons     = sum(idx);
    
    colors_V = getPyPlot_cMap('Dark2',nSesNeurons); colors_V(:,1) = colors_V(:,1).^1.5; colors_V(:,2) = colors_V(:,2).^1.5;
    colors_A = getPyPlot_cMap('Dark2',nSesNeurons); colors_A(:,2) = colors_A(:,2).^1.5; colors_A(:,3) = colors_A(:,3).^1.5;
    
    colors_V = flipud(colors_V);
    colors_A = flipud(colors_A);
    
    tempori         = respmat_ori(idx,:);
    tempfreq        = respmat_freq(idx,:);
    plotlims        = round([0 max([tempfreq(:); tempori(:)])]);
%     plotlims = [0 1];
    figure; set(gcf,'units','normalized','Position',[0.01+(iSes-1)*0.15 0.43 0.3 0.23],'color','w'); hold all
   
    subplot(1,2,1); hold all
    for iN = 1:nSesNeurons
        %         plot(0:45:315,tempori(iN,:),'b-','LineWidth',0.5,'Color',colors_V(iN,:))
        plot(0:45:360,[tempori(iN,:) tempori(iN,1)],'b-','LineWidth',0.5,'Color',colors_V(iN,:))
    end
%     h = shadedErrorBar(0:45:360,[nanmean(tempori(:,:),1) nanmean(tempori(:,1))],[nanstd(tempori(:,:),[],1) nanstd(tempori(:,1))],{'b-','LineWidth',0.5,'Color','k'});
%     delete(h.edge(:));
    
    ylim(plotlims)
    set(gca,'XTick',0:45:360,'YTick',round(get(gca,'ylim')),'FontSize',10)
    xlabel('Orientation'); ylabel('Zscored firing rate response');
    
    subplot(1,2,2); hold all
    for iN = 1:nSesNeurons
        %         plot(8:1:15,tempfreq(iN,:),'r-','LineWidth',0.5,'Color',colors_A(iN,:))
        plot(8:1:16,[tempfreq(iN,:) tempfreq(iN,1)],'r-','LineWidth',0.5,'Color',colors_A(iN,:))
%         plot(8:1:16,[tempfreq(iN,:) tempfreq(iN,1)] / max(tempfreq(iN,:)),'r-','LineWidth',0.5,'Color',colors_A(iN,:))
    end
%     h = shadedErrorBar(8:1:16,[nanmean(tempfreq(:,:),1) nanmean(tempfreq(:,1))],[nanstd(tempfreq(:,:),[],1) nanstd(tempfreq(:,1))],{'b-','LineWidth',0.5,'Color','k'});
%     delete(h.edge(:));

    ylim(plotlims)
    set(gca,'XTick',8:1:16,'YTick',plotlims,'FontSize',10)
    xlabel('Frequency (kHz)'); ylabel('Zscored firing rate response');

end
export_fig(fullfile(params.savedir,'Co-tuning_ExSes_Det'),'-eps','-nocrop')

%% Compute signal correlations between neurons

r_sig_ori   = NaN(nNeurons,nNeurons);
r_sig_frq   = NaN(nNeurons,nNeurons);

for iSes = 1:nSessions
    idx             = find(strcmp(spikeData.session_ID,sessionData.session_ID(iSes)));
    nSesNeurons     = numel(idx);
    
    for iX = 1:nSesNeurons
        for iY = 1:nSesNeurons
            r_sig_ori(idx(iX),idx(iY)) = corr(respmat_ori(idx(iX),:)',respmat_ori(idx(iY),:)');
            r_sig_frq(idx(iX),idx(iY)) = corr(respmat_freq(idx(iX),:)',respmat_freq(idx(iY),:)');
        end
    end
end

r_sig_ori(r_sig_ori>0.99) = NaN; %set auto correlation to NaN
r_sig_frq(r_sig_frq>0.99) = NaN;

%% Make figure of the noise correlation values of tuned neurons:
% this figure can differ from the reported figure and stats because of 
% randomization variability in the shuffle procedure

%Construct indices:
idx_ori             = gOSI_ori_sign==1;
r_sig_ori_tuned     = r_sig_ori(idx_ori,idx_ori);
idx_freq            = gOSI_freq_sign==1;
r_sig_freq_tuned    = r_sig_frq(idx_freq,idx_freq);

%Make figure:
figure; set(gcf,'units','normalized','Position',[0.15 0.43 0.1 0.21],'color','w'); hold all
% bar(1,nanmean(r_sig_ori_tuned(:)),0.6,'b')
% errorbar(1,nanmean(r_sig_ori_tuned(:)),nanstd(r_sig_ori_tuned(:))/sqrt(sum(idx_ori)),'k-','LineWidth',1)
% bar(2,nanmean(r_sig_freq_tuned(:)),0.6,'r')
% errorbar(2,nanmean(r_sig_freq_tuned(:)),nanstd(r_sig_freq_tuned(:))/sqrt(sum(idx_freq)),'k-','LineWidth',1)

datamat  = NaN(500,2);
temp = r_sig_ori_tuned(~isnan(r_sig_ori_tuned));
datamat(1:length(temp),1)= temp;
temp = r_sig_freq_tuned(~isnan(r_sig_freq_tuned));
datamat(1:length(temp),2)= temp;

h   = violinplot(datamat,{'Ori' 'Freq'},'EdgeColor',[1 1 1],'BandWidth',0.1,...
        'Width',0.3,'ShowData',false,'ViolinAlpha',1);
h(1).ViolinColor = 'b';
h(2).ViolinColor = 'r';

xlim([0.5 2.5])
ylim([-0.8 1])
set(gca,'XTick',[1 2],'XTickLabels',{'Ori' 'Freq'},'YTick',[-0.5 0 0.5 1],'FontSize',8)
ylabel('Signal correlation','FontSize',9)

%statistics:
G_mou           = cell(nNeurons,1);
uMice           = unique(sessionData.mousename);
for iMouse = 1:length(uMice)
    G_mou(ismember(spikeData.session_ID,sessionData.session_ID(strcmp(sessionData.mousename,uMice{iMouse})))) = uMice(iMouse);
end

G_mou_ori   = G_mou(gOSI_ori_sign==1);
G_mou_ori   = repmat(G_mou_ori,1,numel(G_mou_ori));
G_mou_freq  = G_mou(gOSI_freq_sign==1);
G_mou_freq  = repmat(G_mou_freq,1,numel(G_mou_freq));

Y_r_sig = [r_sig_ori_tuned(:); r_sig_freq_tuned(:)];
X_mod = [zeros(length(r_sig_ori_tuned(:)),1); ones(length(r_sig_freq_tuned(:)),1)];
G_mou = [G_mou_ori(:); G_mou_freq(:)];

tbl             = table(Y_r_sig,X_mod,G_mou,'VariableNames',{'Rsig','Modality','Mouse'}); %Create table for mixed model
lme             = fitlme(tbl,'Rsig~Modality+(1|Mouse)'); %construct linear mixed effects model with fixed effect of temporal window (NO random intercept for mice as cohorts are tested)
stats           = dataset2table(anova(lme,'DFMethod','Satterthwaite')); %Perform ANOVA on model and output as matrix
fprintf('F(%d,%2.0f)=%1.2f, p=%1.4f; \n',stats{2,3},stats{2,4},stats{2,2},stats{2,5})
p = stats{2,5};
sigstar([1 2],p)
fprintf('n=%d signal correlations from %d orientation-tuned V1 neurons, \n',sum(~isnan(r_sig_ori_tuned(:))),sum(idx_ori))
fprintf('n=%d from %d frequency-tuned V1 neurons)\n',sum(~isnan(r_sig_freq_tuned(:))),sum(idx_freq))

export_fig(fullfile(params.savedir,'SigCor_Violin_Det'),'-eps','-nocrop')
writetable(tbl,'SourceData_S4g_SignalCorr.xlsx')

%%


