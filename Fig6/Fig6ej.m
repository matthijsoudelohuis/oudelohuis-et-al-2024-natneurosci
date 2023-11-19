%% Oude Lohuis et al. 2023 Nat Neurosci
% Triple dissociation of auditory, visual, and motor processing in primary visual cortex
% MOL (C) 2023

startover

%% Parameter settings:
params.Experiments          = {'ChangeDetectionConflictDecor' 'ChangeDetectionConflict'}; %Which versions of the task to load data from
params.ExperimentLabels     = {'NE' 'MST'};         %Labels for the different experiments
params.AlignOn              = 'stimChange';         %On which timestamp to align trials to
params.nExperiments         = length(params.Experiments);
params                      = MOL_getColors_CHDET(params);

params.savedir              = 'E:\Documents\PhD\Figures\Project CHDET\Results - auV1\9ConflictNeural\Correlation';

params                      = params_histresponse_auV1(params);

%Parameters for AUC:
params.min_ntrial           = 10;               %Number of minimum trials that both conditions need to have to be discriminated
params.smoothing            = 0;
params.zscore               = 0;
params.subtr_baseline       = 1;                %Subtract baseline or not
params.nshuffle             = 1000;             %Number of shuffles to base permutation test on
params.alpha                = 0.025;            %Significance level for permutation test

params.minnNeurons          = 10;               %Minimum number of simultaneously recorded neurons to compute fraction of tuned cells


%Set colors and labels for conditions:
params.colors_splits = [];
params.colors_splits(1,1:2,:) = [160 160 160; 143 43 44];% 157 0 214]; 
params.colors_splits(2,1:2,:) = [37 45 138; 157 0 214];% 157 0 214]; 
params.colors_splits = params.colors_splits / 256;
params.labels_splits = {'av' 'Av'; 'aV' 'AV'};
params.lines_splits = {'-'  '--'; '-.' ':'};

%% Get data:
[Data] = MOL_GetData('E:','CHDET',params.Experiments,{},[],{'sessionData' 'trialData_newtrials' 'spikeData'});
sessionData         = Data.sessionData;
trialData           = Data.trialData;
spikeData           = Data.spikeData;

%% Remove last n trials:
trialData           = MOL_RemoveLastnTrials(trialData,5);

%% Filter out neurons based on quality:
spikeData               = MOL_filterNeurons(sessionData,trialData,spikeData);

%% Filter out sessions with muscimol:
sesids              = sessionData.session_ID(cellfun(@isempty,sessionData.MuscimolArea));
fprintf('Removed %d/%d sessions  with muscimol manipulations\n',length(sessionData.session_ID)-length(sesids),length(sessionData.session_ID));
[sessionData, trialData, spikeData]        = MOL_getTempPerSes(sesids,sessionData,trialData,spikeData);

%% Filter out trials with photostimulation in V1:
sesids              = sessionData.session_ID(sessionData.UseOpto & strcmp(sessionData.PhotostimArea,'V1'));
idx                 = trialData.hasphotostim==1 & ismember(trialData.session_ID,sesids);
trialFields         = fieldnames(trialData);
fprintf('Removed %d/%d trials  with photostim in V1\n',sum(idx),length(trialData.session_ID));
for iF = 1:length(trialFields)
    trialData.(trialFields{iF}) = trialData.(trialFields{iF})(~idx);
end

%% Filter neurons in area:
idx             = ismember(spikeData.area,{'V1'});
spikeFields     = fieldnames(spikeData);
fprintf('Subselected %d V1 neurons \n',sum(idx));
for iF = 1:length(spikeFields)
    spikeData.(spikeFields{iF}) = spikeData.(spikeFields{iF})(idx);
end

%% Filter sessions that have one or four postFreq and postOri during conflict trials:
nSessions           = length(sessionData.session_ID);
singleStimConfl     = false(nSessions,1);
for iSes = 1:nSessions
    postfreqs               = trialData.audioFreqPostChangeNorm(strcmp(trialData.trialType,'C') & strcmp(trialData.session_ID,sessionData.session_ID(iSes)));
    postoris                = trialData.visualOriPostChangeNorm(strcmp(trialData.trialType,'C') & strcmp(trialData.session_ID,sessionData.session_ID(iSes)));
    singleStimConfl(iSes)   = ismember(numel(unique(postfreqs)),[1 4]) && ismember(numel(unique(postoris)),[1 4]);
end
fprintf('Removed %d/%d sessions with multiple conflict stimuli\n',length(sessionData.session_ID)-sum(singleStimConfl),length(sessionData.session_ID));
[sessionData,trialData,spikeData] = MOL_getTempPerSes(sessionData.session_ID(singleStimConfl),sessionData,trialData,spikeData);

%% Get only sessions with neurons recorded:
[sessionData,trialData,spikeData] = MOL_getTempPerSes(unique(spikeData.session_ID),sessionData,trialData,spikeData);

%% Save dataset:
save('Dataset6_1.mat','params','sessionData','trialData','spikeData')

%% Or load dataset
load Dataset6_1.mat

%% Report dataset:
fprintf('Dataset: %d sessions, %d trials, %d neurons\n',length(sessionData.session_ID),length(trialData.session_ID),length(spikeData.session_ID));







%% Compute Z-scored psth for all combinations of modalities and saliencies:

%Parameters for the firing rate calculations:
params                  = params_histresponse_auV1(params);

nNeurons                = length(spikeData.ts); %number of neurons in total (both cohorts, all animals, all session)
nVsals                  = 3; %number of visual saliencies (thr and max)
nAsals                  = 3; %number of auditory saliencies (thr and max)

nTrialtypes             = 3;
zmat                    = NaN(nNeurons,params.nTimebins,nVsals,nAsals); %init variables (Rvs is correlation of firing rate during conflict trial to stimulus-matched visual unimodal trial)

params.conv_win         = 'chg';

fprintf('Computing firing rate for neuron        \n');

for iNeuron = 1:nNeurons %Loop over all neurons:
    fprintf(repmat('\b', 1, numel([num2str(iNeuron-1) num2str(nNeurons)])+2));
    fprintf('%d/%d\n',iNeuron,nNeurons);
    
    %Get the relevant data for each session individually:
    temptrialData        = MOL_getTempPerSes(spikeData.session_ID(iNeuron),trialData);
    
    %Compute histogram:
    events_ts                               = temptrialData.(params.AlignOn);
    hist_mat                                = calc_psth(events_ts,spikeData.ts{iNeuron},params);    %Construct histogram matrix:
    
    temptrialData.visualOriPostChangeNorm   = round(temptrialData.visualOriPostChangeNorm/2);
    temptrialData.audioFreqPostChangeNorm   = round(temptrialData.audioFreqPostChangeNorm/2);
    
    for iVsal = 1:nVsals
        for iAsal = 1:nAsals
            idx   = temptrialData.visualOriChangeNorm==iVsal & temptrialData.audioFreqChangeNorm==iAsal;
            
            zmat(iNeuron,:,iVsal,iAsal)                = nanmean(hist_mat(idx,:),1);
        end
    end
end

%%
% idx_time = params.xtime>0 & params.xtime<500e3;
% 
% zmat2 = reshape(zmat,nNeurons,params.nTimebins,nTrialtypes,nVsals*nAsals);
% 
% tempresp = squeeze(nanmean(nanmean(zmat2(:,idx_time,:,:,:),2),4));
% 
% [~,sortidx] = sort(tempresp(:,1) - tempresp(:,2));
% 
% zmat2 = zmat2(sortidx,:,:,:);

%% Fig 6e - 3x3 heatmap of z-scored activity across the population

params.labels_splits = {'C' 'a-' 'A-'; 'v-' 'av' 'Av'; 'V-' 'aV' 'AV'};

for iExp = [1 2]
    figure; hold all; set(gcf,'units','normalized','Position',[0.15 0.15 0.3 0.5],'color','w');
    idx_exp         = ismember(spikeData.session_ID,sessionData.session_ID(strcmp(sessionData.Experiment,params.Experiments(iExp))));
    
%      idx_exp         = ismember(spikeData.session_ID,sessionData.session_ID(strcmp(sessionData.Experiment,params.Experiments(iExp)))) & ...
%         ~(ismember(spikeData.session_ID,sessionData.session_ID(ismember(sessionData.mousename,{'2044' '2045'})))); %Get all neurons from specific animals (NE, UST, MST)
    
    %Sort the zmat:
    idx_time        = params.xtime>0e3 & params.xtime<500e3;
    zmat2           = zmat(idx_exp,:,:,:);
    tempVresp       = squeeze(nanmean(nanmean(zmat2(:,idx_time,2:3,1),3),2));
    tempAresp       = squeeze(nanmean(nanmean(zmat2(:,idx_time,1,2:3),4),2));
    [~,sortidx]     = sort(tempVresp - tempAresp);
%     [~,sortidx]     = sort(tempAresp - tempVresp);
    zmat2           = zmat2(sortidx,:,:,:,:);
    
    for iVsal = 1:nVsals
        for iAsal = 1:nAsals
%             subplot(nVsals,nAsals,(iVsal-1)*nAsals + iAsal); hold all;
                        subplot(nVsals,nAsals,(iAsal-1)*nVsals + iVsal); hold all;
            tempdat = squeeze(zmat2(:,:,iVsal,iAsal));
            nNeuronsExp = sum(idx_exp);
            imagesc(params.xtime,1:nNeuronsExp,tempdat);
            plot([0 0],[0 nNeuronsExp],':','Color',[1 1 1],'LineWidth',1)
            caxis([-0.5 3])
            title(params.labels_splits{iVsal,iAsal})
            xlim([params.t_pre params.t_post])
            ylim([1 nNeuronsExp])
            if iAsal==3
                set(gca,'XTick',(-4:0.5:4)*1e6,'XTickLabels',(-4:0.5:4)*1e3);
                xlabel('Time (ms)')
            else
                set(gca,'XTick',[])
            end
            if iVsal==1 && iAsal==2
                set(gca,'YTick',[1 nNeuronsExp]); ylabel('Neurons')
            else
                set(gca,'YTick',[])
            end
        end
    end
    tightfig()
    filename = sprintf('Zmat_3x3_%s.eps',params.ExperimentLabels{iExp});
    export_fig(fullfile(params.savedir,filename),gcf);
end



%% AUC Analysis

%Parameters for AUC:
params.min_ntrial           = 10;               %Number of minimum trials that both conditions need to have to be discriminated
params.smoothing            = 0;
params.zscore               = 0;
params.subtr_baseline       = 1;                %Subtract baseline or not
params.nshuffle             = 1;             %Number of shuffles to base permutation test on
params.alpha                = 0.025;            %Significance level for permutation test

params.twin_resp_start      = 0e6;
params.twin_resp_stop       = 0.5e6;

%% Compute AUC for neurons:
[AUC_ORI_V,AUC_ORI_C,pVal_ORI_V,pVal_ORI_C,respVis,respAud] = calc_AUC_spikes_conflict(params,sessionData,trialData,spikeData);

%% Fig 6j - correlation between AUC during visual and audiovisual trials:

iExp            = 2;
idx_exp         = ismember(spikeData.session_ID,sessionData.session_ID(strcmp(sessionData.Experiment,params.Experiments(iExp))));

idx_V1          = strcmp(spikeData.area,'V1') & idx_exp;


figure; set(gcf,'units','normalized','Position',[0.25 0.37 0.23 0.39],'color','w'); hold all
sigori = pVal_ORI_V<0.025 | pVal_ORI_C<0.025 | pVal_ORI_V>0.975 | pVal_ORI_C>0.975;
scatter(AUC_ORI_V(sigori & idx_V1),AUC_ORI_C(sigori & idx_V1),25,[0 0 0.8],'filled');
scatter(AUC_ORI_V(~sigori & idx_V1),AUC_ORI_C(~sigori & idx_V1),25,[0.4 0.4 0.4],'filled');
xlim([-1 1]); ylim([-1 1]);
plot([0 0],[-1 1],'k:','LineWidth',0.5)
plot([-1 1],[0 0],'k:','LineWidth',0.5)
set(gca,'XTick',[-1 0 1],'YTick',[-1 0 1]);
xlabel('Visual'); ylabel('Audiovisual')
title('Orientation')

nNeurons = length(spikeData.session_ID);
G_mou           = cell(nNeurons,1);
uMice           = unique(sessionData.mousename);
for iMouse = 1:length(uMice)
    G_mou(ismember(spikeData.session_ID,sessionData.session_ID(strcmp(sessionData.mousename,uMice{iMouse})))) = uMice(iMouse);
end

tbl             = table(AUC_ORI_C(idx_V1),AUC_ORI_V(idx_V1),G_mou(idx_V1),'VariableNames',{'AUC_ORI_C','AUC_ORI_V','Mouse'}); %Create table for mixed model
lme             = fitlme(tbl,'AUC_ORI_C~AUC_ORI_V+(1|Mouse)'); %construct linear mixed effects model with fixed effect of temporal window and random intercept for different mice
stats           = dataset2table(anova(lme,'DFMethod','Satterthwaite')); %Perform ANOVA on model and output as matrix
p = stats{2,5};
fprintf('AUC are strongly correlated (R=%1.3f, F(%d,%2.0f)=%1.2f, p=%1.2e; \n',sqrt(lme.Rsquared.Ordinary),stats{2,3},stats{2,4},stats{2,2},stats{2,5})
if p<0.05
    sigstar([0.5 1],p)
end

MOL_prepfigAI

export_fig(fullfile(params.savedir,sprintf('Scatter_AUC_V_AV')),'-eps','-nocrop')

writetable(tbl,'SourceData_Fig6j_AUC_V_AV.xlsx')



%%









