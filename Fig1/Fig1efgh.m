%% This script describes some elementary features of auditory responses in V1 during 
% an audiovisual change detection task for three different task contingencies 
% as reported in Oude Lohuis et al. 2024

%% Parameter settings:

params.Experiments          = {'ChangeDetectionConflictDecor' 'VisOnlyTwolevels' 'ChangeDetectionConflict'}; %Which versions of the task to load data from
params.ExperimentLabels     = {'NE' 'UST' 'MST'}; %Which versions of the task to load data from
params.AlignOn              = 'stimChange';         %On which timestamp to align trials to
params.nExperiments         = length(params.Experiments);
params                      = MOL_getColors_CHDET(params);

params                      = params_histresponse_auV1(params);

% Statistics:
params.alpha                = 0.025;
params.posthoctest          = 'bonferroni';

params.areas                = {'V1'};
params.nAreas               = length(params.areas);

params.minTrialCond         = 10; %number of trial that need to be present to calculate mean response and test for significant response
params.minNfrac             = 15; %number of simultaneously recorded neurons to compute fraction responsive

params.nshuffle             = 1000; %Number of shuffles to base permutation test on
params.alpha                = 0.05; %Significance level for permutation test

params.savedir              = 'E:\Documents\PhD\Figures\Project CHDET\Results - auV1\2AuV1recruitment';

%% Get data:
[Data] = MOL_GetData('E:','CHDET',params.Experiments,{},[],{'sessionData' 'trialData_newtrials' 'spikeData'});
sessionData         = Data.sessionData;
trialData           = Data.trialData;
spikeData           = Data.spikeData;

%% Remove last 20 trials:
trialData           = MOL_RemoveLastnTrials(trialData,20);

%% Filter out neurons based on quality:
spikeData           = MOL_filterNeurons(sessionData,trialData,spikeData);

%% Remove sessions that are passive/active from 2044 and 2045s: 
sesids                      = sessionData.session_ID(~ismember(sessionData.mousename,{'2044' '2045'}));
fprintf('Removed %d/%d sessions from 2044 and 2045 with active and passive epochs\n',length(sessionData.session_ID)-length(sesids),length(sessionData.session_ID));
[sessionData, trialData, spikeData]        = MOL_getTempPerSes(sesids,sessionData,trialData,spikeData);

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

%% Filter neurons on area:
spikeFields     = fieldnames(spikeData);
idx             = ismember(spikeData.area,params.areas);
fprintf('Filtered %d/%d neurons based on area\n',sum(idx),length(spikeData.session_ID));
for iF = 1:length(spikeFields)
    spikeData.(spikeFields{iF}) = spikeData.(spikeFields{iF})(idx);
end

%% Get only sessions with neurons recorded:
sesids              = unique(spikeData.session_ID);
fprintf('Removed %d/%d sessions without neurons in target area\n',length(sessionData.session_ID)-length(sesids),length(sessionData.session_ID));
[sessionData,trialData,spikeData] = MOL_getTempPerSes(sesids,sessionData,trialData,spikeData);

%% Save dataset:
save('Dataset1_2.mat','params','sessionData','trialData','spikeData')

%% Or start script from saved dataset:
load Dataset1_2.mat

%% Report dataset:
fprintf('Dataset: %d sessions, %d trials, %d neurons\n',length(sessionData.session_ID),length(trialData.session_ID),length(spikeData.session_ID));

%%




%% Fig 1e Plot some rasters for example cells:
%Auditory responsive example neurons:
cell_IDs = {'10092019030831134'     %NE
    '20292019121111234'             %UST
    '20122018081431146'};           %MST

params.labels_neurons       = params.ExperimentLabels;

params.t_pre                = -0.25e6;
params.t_post               = 0.75e6;

%Construct bin edges and time axis
params.edges                = [params.t_pre:params.binsize:params.t_post] - params.binsize/2;                    %#ok<NBRAK> Define edges of the bins
params.xtime                = params.t_pre:params.binsize:params.t_post-params.binsize;                        %Define time axis
params.nTimebins            = length(params.xtime); %number of time bins

params.zscore               = 0;

params.exportfig            = 0;

plotRaster_AuResp(sessionData,trialData,spikeData,cell_IDs,params)

%% Init output fields:
nSessions                       = length(sessionData.session_ID);
nNeurons                        = length(spikeData.session_ID);

spikeData.sign_visresponse      = NaN(nNeurons,1);
spikeData.sign_auresponse       = NaN(nNeurons,1);

spikeData.z_visresponse         = NaN(nNeurons,1);
spikeData.z_auresponse          = NaN(nNeurons,1);

%% For each neuron compute sign resp to idx:
params                          = params_histresponse_auV1(params); %get standard parameters for AuV1 paper again

params.nSplits                  = 4; %four trial types (2 per modality)

%Longer baseline to compute stable baseline activity:
params.t_pre                    = -1e6;
params.edges                    = [params.t_pre:params.binsize:params.t_post] - params.binsize/2;                    %#ok<NBRAK> Define edges of the bins
params.xtime                    = params.t_pre:params.binsize:params.t_post-params.binsize;                        %Define time axis
params.nTimebins                = length(params.xtime); %number of time bins
twin_baseline_start             = -1e6;
twin_baseline_stop              = 0;

lastsesid                       = []; %memory var to keep track of whether neuron comes from different session and design matrix needs to be reconstructed
zmat                            = NaN(nNeurons,params.nTimebins,params.nSplits); %init output var

fprintf('Computing average Z-scored response for neuron        \n');

for iNeuron = 1:nNeurons %Loop over all neurons:
    fprintf(repmat('\b', 1, numel([num2str(iNeuron-1) num2str(nNeurons)])+2));
    fprintf('%d/%d\n',iNeuron,nNeurons);
    
    if ~strcmp(lastsesid,spikeData.session_ID(iNeuron)) %construct new predictor matrix if neuron comes from a new session:
        %Get the relevant data for each session individually:
        temptrialData        = MOL_getTempPerSes(spikeData.session_ID(iNeuron),trialData);
        lastsesid            = spikeData.session_ID(iNeuron); %save this session_ID
    end
    
    %Compute histogram:
    events_ts               = temptrialData.stimChange;
    hist_mat                = calc_psth(events_ts,spikeData.ts{iNeuron},params);    %Construct histogram matrix:
    
    %identify the splits
    splits                  = {};
    splits{1}               = ismember(temptrialData.trialType,{'X'}) & temptrialData.visualOriChangeNorm==3 & ismember(temptrialData.visualOriPostChangeNorm,[1 2]);
    splits{2}               = ismember(temptrialData.trialType,{'X'}) & temptrialData.visualOriChangeNorm==3 & ismember(temptrialData.visualOriPostChangeNorm,[3 4]);
    splits{3}               = ismember(temptrialData.trialType,{'Y'}) & temptrialData.audioFreqChangeNorm==3 & ismember(temptrialData.audioFreqPostChangeNorm,[1 2]);
    splits{4}               = ismember(temptrialData.trialType,{'Y'}) & temptrialData.audioFreqChangeNorm==3 & ismember(temptrialData.audioFreqPostChangeNorm,[3 4]);
    
    bsl                     = nanmean(hist_mat(:,params.xtime>params.twin_baseline_start & params.xtime<params.twin_baseline_stop),2);
    resp                    = nanmean(hist_mat(:,params.xtime>params.twin_resp_start & params.xtime<params.twin_resp_stop),2);
    
    splitsign               = false(4,1);
    for iSplit = 1:params.nSplits %Store the mean response for each of these splits
        if sum(splits{iSplit})>=params.minTrialCond && ~all(isnan(bsl))
            zmat(iNeuron,:,iSplit)  = nanmean(hist_mat(splits{iSplit},:),1);
            
            [~,splitsign(iSplit)]   = signrank(bsl(splits{iSplit}),resp(splits{iSplit}),'alpha',params.alpha);
        end
    end
    
    spikeData.sign_visresponse(iNeuron,1)       = any(splitsign(1:2)); %significant if either of the stimul is significant.
    spikeData.sign_auresponse(iNeuron,1)        = any(splitsign(3:4));
    
    spikeData.z_visresponse(iNeuron,1)          = max([nanmean(resp(splits{1})) nanmean(resp(splits{2}))]);
    spikeData.z_auresponse(iNeuron,1)           = max([nanmean(resp(splits{3})) nanmean(resp(splits{4}))]);
    
    %incorporate possibility for reduction of firing rate:
    if abs(min([nanmean(resp(splits{1})) nanmean(resp(splits{2}))])) > abs(max([nanmean(resp(splits{1})) nanmean(resp(splits{2}))]))
        spikeData.z_visresponse(iNeuron,1)          = min([nanmean(resp(splits{1})) nanmean(resp(splits{2}))]);
    else
        spikeData.z_visresponse(iNeuron,1)          = max([nanmean(resp(splits{1})) nanmean(resp(splits{2}))]);
    end
    if abs(min([nanmean(resp(splits{3})) nanmean(resp(splits{4}))])) > abs(max([nanmean(resp(splits{3})) nanmean(resp(splits{4}))]))
        spikeData.z_auresponse(iNeuron,1)          = min([nanmean(resp(splits{3})) nanmean(resp(splits{4}))]);
    else
        spikeData.z_auresponse(iNeuron,1)          = max([nanmean(resp(splits{3})) nanmean(resp(splits{4}))]);
    end
    
end

%% Figure 1f that shows the mean z-scored response over all positively responding auditory V1 neurons

%Get mean response during 0-200 ms:
zresp = squeeze(nanmean(zmat(:,params.xtime>params.twin_resp_start & params.xtime<=params.twin_resp_stop,:),2));

%Get the response for the preferred auditory stimulus:
[~,idx_max]                 = max(zresp(:,[3 4]),[],2);
matAudRate                  = NaN(nNeurons,params.nTimebins);
matAudRate(idx_max==1,:)    = zmat(idx_max==1,:,3);
matAudRate(idx_max==2,:)    = zmat(idx_max==2,:,4);

%Make figure:
figure; hold all; set(gcf,'units','normalized','Position',[0.12 0.37 0.16 0.2],'color','w');
handles = [];
fprintf('n=')

for iExp = 1:params.nExperiments %for each experiment select appropriate neurons idx
    idx_exp         = ismember(spikeData.session_ID,sessionData.session_ID(strcmp(sessionData.Experiment,params.Experiments(iExp))));
    idx_resp        = spikeData.sign_auresponse==1;% & spikeData.z_auresponse>0;
    idx_all         = idx_exp & idx_resp;
    fprintf('%d, ',sum(idx_all))
    handles(iExp)   = plot(params.xtime*1e-3,squeeze(nanmean(matAudRate(idx_all,:),1)),'Color',params.colors_experiments{iExp}); %#ok<SAGROW>
end
fprintf(' auditory responsive neurons for NE, UST, MST respectively\n')
% rectangle('Position',[0 0 200 2.25],'FaceColor',[0.3 0.3 .3],'EdgeColor','none'); %uncomment for shaded window display

%Figure make up:
xlim([-0.25 0.75]*1e3)
ylim([-0.25 3])
plot([0 0],[-0.5 1.75],'-','Color',[0.6 0.6 0.6],'LineWidth',0.25)
legend(handles,params.ExperimentLabels,'Location','NorthEast'); legend boxoff;
set(gca,'XTick',[-250 0 250 500 750],'YTick',[0 1 2 3])
xlabel('Time (ms)')

MOL_prepfigAI
filename = sprintf('V1_auposresp_ztime_mean_3cohorts.eps');
% export_fig(fullfile(params.savedir,filename),gcf);

%% Figure 1f inset: Mean z-scored response of auditory responsive neurons in V1 across the cohorts:
nInit               = 1000;
datatoplot          = NaN(params.nExperiments,nInit);
idx_resp            = spikeData.sign_auresponse==1;% & spikeData.z_auresponse>0;

for iExp = 1:params.nExperiments
    idx_exp         = ismember(spikeData.session_ID,sessionData.session_ID(strcmp(sessionData.Experiment,params.Experiments(iExp))));
    idx_all         = idx_exp & idx_resp;
    
    datatoplot(iExp,1:sum(idx_all)) = spikeData.z_auresponse(idx_all);
    fprintf('% \n')
end

%Compute mean and sem:
meantoplot      = nanmean(datatoplot,2);
errortoplot     = nanstd(datatoplot,[],2) ./ sum(~isnan(datatoplot),2);

%Make the figure;
figure; hold all; set(gcf,'units','normalized','Position',[0.4 0.58 0.06 0.19],'color','w');

barlocations = [1 1.5 2];
z = errorbar(barlocations, meantoplot, errortoplot,'r.','MarkerSize',15,'LineWidth',2);
errorbar_tick(z,0.001,'units')

%Multi-level statistics: 
X_coh           = cell(nNeurons,1);
X_coh(ismember(spikeData.session_ID,sessionData.session_ID(ismember(sessionData.Experiment,params.Experiments{1}))))         = params.ExperimentLabels(1);
X_coh(ismember(spikeData.session_ID,sessionData.session_ID(ismember(sessionData.Experiment,params.Experiments{2}))))         = params.ExperimentLabels(2);
X_coh(ismember(spikeData.session_ID,sessionData.session_ID(ismember(sessionData.Experiment,params.Experiments{3}))))         = params.ExperimentLabels(3);

Y               = spikeData.z_auresponse;

G_mou           = cell(nNeurons,1);
uMice           = unique(sessionData.mousename);
for iMouse = 1:length(uMice)
    G_mou(ismember(spikeData.session_ID,sessionData.session_ID(strcmp(sessionData.mousename,uMice{iMouse})))) = uMice(iMouse);
end

tbl             = table(Y(idx_resp),X_coh(idx_resp),G_mou(idx_resp),'VariableNames',{'Activity','Cohort','Mouse'}); %Create table for mixed model

lme             = fitlme(tbl,'Activity~Cohort'); %construct linear mixed effects model with fixed effect of temporal window (NO random intercept for mice as cohorts are tested)
stats           = dataset2table(anova(lme,'DFMethod','Satterthwaite')); %Perform ANOVA on model and output as matrix
fprintf('Z-scored auditory activity different between cohorts: \n')
fprintf('F(%d,%2.0f)=%1.2f, p=%1.2f; \n',stats{2,3},stats{2,4},stats{2,2},stats{2,5})

contrasts       = {[0 1 0] [0 -1 1] [0 0 1]}; %which contrasts to do (coding is relative to first category if only one dummy)
contrastpos     = {[2 1.5] [1.5 2] [1 2]};
contrastlabels = {'MST vs NE' 'NE vs UST' 'MST vs UST'};

fprintf('Posthoc comparison:\n')
for iC = 1:3
    [p,F,DF1,DF2] = coefTest(lme,contrasts{iC});
    fprintf('(%s):F(%d,%d)=%1.2f, p=%1.3f \n',contrastlabels{iC},DF1,DF2,F,p)
    if p<0.05
        sigstar(contrastpos{iC},p);
    end
end

%Figure make up:
xlim([0.7 2.3])
ylim([-.5 2.5])
set(gca,'XTick',barlocations,'XTickLabel',params.ExperimentLabels,'XTickLabelRotation',45,'Fontsize',8);
ylabel('Auditory response (z-score)','FontSize',10)
%Export figure:
filename = sprintf('MeanZAu_V1_bar_3cohorts.eps');
% export_fig(fullfile(params.savedir,filename),gcf);
writetable(tbl,'SourceData_Fig1f_auresp_cohorts.xlsx')

%% Convert significance value to logical:
spikeData.sign_visresponse      = spikeData.sign_visresponse==1;
spikeData.sign_auresponse       = spikeData.sign_auresponse==1;

%Compute for each session the number of auditory and visually responsive units in V1:
FracVisResp_V1         = NaN(nSessions,1);
FracAudResp_V1         = NaN(nSessions,1);
MeanZVisResp_V1        = NaN(nSessions,1);
MeanZAudResp_V1        = NaN(nSessions,1);

for iSes = 1:nSessions %loop over sessions
    idx_ses             = strcmp(spikeData.session_ID,sessionData.session_ID(iSes));
    if sum(idx_ses)>=params.minNfrac
        FracVisResp_V1(iSes)   = sum(spikeData.sign_visresponse(idx_ses))/sum(idx_ses);
        FracAudResp_V1(iSes)   = sum(spikeData.sign_auresponse(idx_ses))/sum(idx_ses);
        %Store mean response for this session as well:
        MeanZVisResp_V1(iSes)  = nanmean(spikeData.z_visresponse(idx_ses & spikeData.sign_visresponse));
        MeanZAudResp_V1(iSes)  = nanmean(spikeData.z_auresponse(idx_ses & spikeData.sign_auresponse));
    end
end

%% Show fraction of auditory responsive neurons in V1 across the cohorts:
nInit               = 1000;
datatoplot          = NaN(params.nExperiments,nInit);

for iExp = 1:params.nExperiments
    datatoplot(iExp,1:sum(strcmp(sessionData.Experiment,params.Experiments(iExp)))) = FracAudResp_V1(strcmp(sessionData.Experiment,params.Experiments(iExp)));
%     datatoplot(iExp,1:sum(strcmp(sessionData.Experiment,params.Experiments(iExp)))) = FracVisResp(strcmp(sessionData.Experiment,params.Experiments(iExp)));
end

%Compute mean and sem:
meantoplot      = nanmean(datatoplot,2);
errortoplot     = nanstd(datatoplot,[],2) ./ sum(~isnan(datatoplot),2);

%Make the figure;
figure; hold all; set(gcf,'units','normalized','Position',[0.4 0.58 0.1 0.18],'color','w');

barlocations = 1:3;
z = errorbar(1:3, meantoplot, errortoplot,'r.','MarkerSize',15,'LineWidth',2);
errorbar_tick(z,0.001,'units')

set(gca,'XTick',1:3,'XTickLabel',params.ExperimentLabels,'XTickLabelRotation',45,'Fontsize',8);
ylabel('% Responsive neurons','FontSize',10)

xlim([0.5 3.5])
ylim([0 0.6])

%Statistical testing:
groups              = repmat(1:3,nInit,1); groups = reshape(groups,params.nExperiments*nInit,1);
datatotest          = permute(datatoplot,[2 1]); %get number of inits first
datatotest          = reshape(datatotest,params.nExperiments*nInit,1); %reshape to one column vector
groups              = groups(~isnan(datatotest)); %filter out nans
datatotest          = datatotest(~isnan(datatotest)); %filter out nans
%perform kruskal wallis nonparametric anova:
[p,~,stats]         = kruskalwallis(datatotest,groups,'off');
comptable           = multcompare(stats,'display','off','alpha',0.05,'ctype',params.posthoctest); %do posthoc
comptable(:,[1 2])  = barlocations(comptable(:,[1 2])); %replace groups by bar locations
comptable           = comptable(comptable(:,end)<0.05,:); %Filter only significant
sigstar(mat2cell(comptable(:,1:2),ones(size(comptable,1),1)),comptable(:,end)) %use sigstar function to identify

fprintf('Fraction of auditory responsive neurons in V1 per cohort:\n')
for iExp = 1:params.nExperiments
    fprintf('%s: %2.1f +/- %2.1f %%,n=%d sessions;\n',params.ExperimentLabels{iExp},meantoplot(iExp)*100,errortoplot(iExp)*100,sum(~isnan(datatoplot(iExp,:))))
end
fprintf('p=%1.2f, Kruskal-Wallis test)\n',p)

filename = sprintf('FracAu_V1_bar_3cohorts.eps');
export_fig(fullfile(params.savedir,filename),gcf);

tbl             = table(datatotest,groups,'VariableNames',{'Frac_au_Responsive','Cohorts'}); %Create table for mixed model
writetable(tbl,'SourceData_Frac_Auresposive.xlsx')

%% Show fraction of visually responsive neurons across the cohorts:
nInit               = 1000;
datatoplot          = NaN(params.nExperiments,nInit);

for iExp = 1:params.nExperiments
%     datatoplot(iExp,1:sum(strcmp(sessionData.Experiment,params.Experiments(iExp)))) = FracAudResp(strcmp(sessionData.Experiment,params.Experiments(iExp)));
    datatoplot(iExp,1:sum(strcmp(sessionData.Experiment,params.Experiments(iExp)))) = FracVisResp_V1(strcmp(sessionData.Experiment,params.Experiments(iExp)));
end

%Compute mean and sem:
meantoplot      = nanmean(datatoplot,2);
errortoplot     = nanstd(datatoplot,[],2) ./ sum(~isnan(datatoplot),2);

%Make the figure;
figure; hold all; set(gcf,'units','normalized','Position',[0.4 0.58 0.1 0.18],'color','w');

barlocations = 1:3;
z = errorbar(1:3, meantoplot, errortoplot,'b.','MarkerSize',15,'LineWidth',2);
errorbar_tick(z,0.001,'units')

set(gca,'XTick',1:3,'XTickLabel',params.ExperimentLabels,'XTickLabelRotation',45,'Fontsize',8);
ylabel('% Responsive neurons','FontSize',10)

xlim([0.5 3.5])
ylim([0 0.6])

%Statistical testing:
groups              = repmat(1:3,nInit,1); groups = reshape(groups,params.nExperiments*nInit,1);
datatotest          = permute(datatoplot,[2 1]); %get number of inits first
datatotest          = reshape(datatotest,params.nExperiments*nInit,1); %reshape to one column vector
groups              = groups(~isnan(datatotest)); %filter out nans
datatotest          = datatotest(~isnan(datatotest)); %filter out nans
%perform kruskal wallis nonparametric anova:
[p,~,stats]         = kruskalwallis(datatotest,groups,'off');
comptable           = multcompare(stats,'display','off','alpha',0.05,'ctype',params.posthoctest); %do posthoc
comptable(:,[1 2])  = barlocations(comptable(:,[1 2])); %replace groups by bar locations
comptable           = comptable(comptable(:,end)<0.05,:); %Filter only significant
sigstar(mat2cell(comptable(:,1:2),ones(size(comptable,1),1)),comptable(:,end)) %use sigstar function to identify

disp(comptable(:,[1 2 end]))

fprintf('Fraction of visually responsive neurons in V1 per cohort:\n')
for iExp = 1:params.nExperiments
    fprintf('%s: %2.1f +/- %2.1f %%,n=%d sessions;\n',params.ExperimentLabels{iExp},meantoplot(iExp)*100,errortoplot(iExp)*100,sum(~isnan(datatoplot(iExp,:))))
end
fprintf('p=%1.2f, Kruskal-Wallis test)\n',p)

filename = sprintf('FracVis_V1_bar_3cohorts.eps');
export_fig(fullfile(params.savedir,filename),gcf);

tbl             = table(datatotest,groups,'VariableNames',{'Frac_Vis_Responsive','Cohorts'}); %Create table for mixed model
writetable(tbl,'SourceData_Frac_Visresposive.xlsx')

%%







%% %%%%%%%%%%
%Frequency tuning: 

%% Figure 1h: Plot some rasters for tuned example cells:
%Auditory responsive example neurons:
cell_IDs = {
    '10092019030831126'     %NE
    '20122018081311186'     %MST
    };

params.labels_neurons       = cell_IDs;

params.t_pre                = -0.25e6;
params.t_post               = 0.75e6;

%Construct bin edges and time axis
params.edges                = [params.t_pre:params.binsize:params.t_post] - params.binsize/2;                    %#ok<NBRAK> Define edges of the bins
params.xtime                = params.t_pre:params.binsize:params.t_post-params.binsize;                        %Define time axis
params.nTimebins            = length(params.xtime); %number of time bins

params.zscore               = 0;

params.exportfig            = 1;

params.colors_splits = {[0.2 0.6 0.4] [0.4 0.3 0.9]};

plotRaster_AuFreq(sessionData,trialData,spikeData,cell_IDs,params)

%% Compute AUC for neurons: *takes quite some time
[AUC_ORI_CELL,AUC_FREQ_CELL,pVal_ORI_CELL,pVal_FREQ_CELL,respVis,respAud] = calc_AUC_spikes(params,sessionData,trialData,spikeData);

%% Binarize significantly tuned neurons:
sigOri                      = pVal_ORI_CELL<params.alpha/2 | pVal_ORI_CELL>(1-params.alpha/2); %divide by two for two-tailed test
sigFreq                     = pVal_FREQ_CELL<params.alpha/2 | pVal_FREQ_CELL>(1-params.alpha/2);

%% Quantification of fraction significantly tuned neurons:
fracOriSes              = NaN(nSessions,1);
fracFreqSes             = NaN(nSessions,1);

for iSes = 1:nSessions
    idx = strcmp(spikeData.session_ID,sessionData.session_ID(iSes));% & spikeData.celltype==2;
    if sum(idx)>=params.minNfrac && sum(~isnan(AUC_ORI_CELL(idx)))>=params.minNfrac
        fracOriSes(iSes) = sum(sigOri(idx)) / sum(idx);
        fracFreqSes(iSes) = sum(sigFreq(idx)) / sum(idx);
    end
end

%% Report the fraction of tuned neurons per session:
temp = [nanmean(fracFreqSes) nanstd(fracFreqSes)/sqrt(sum(~isnan(fracFreqSes))) nanmean(fracOriSes) nanstd(fracOriSes)/sqrt(sum(~isnan(fracOriSes)))]*100;
fprintf('Fraction of V1 neurons that sound-evoked firing rate (0-200 ms) was stimulus specific\n was %2.1f +- %1.1f (Mean +- SEM across sessions), which was similar to the fraction\n of orientation selective cells (%2.1f +- %1.1f)\n',temp)

%% Report the fraction of tuned neurons per session:
fprintf('Fraction of frequency tuned neurons in V1 per cohort:\n')
for iExp = 1:params.nExperiments
    idx = strcmp(sessionData.Experiment,params.Experiments{iExp});
    temp = [nanmean(fracFreqSes(idx)) nanstd(fracFreqSes(idx)) / sqrt(sum(idx))]*100;
    fprintf('%s: %2.1f +- %1.1f\n',params.ExperimentLabels{iExp},temp)
end

[p,~,stats]         = kruskalwallis(fracFreqSes,sessionData.Experiment,'off');

fprintf('p=%1.2f, Kruskal-Wallis test)\n',p)

%% Statistical difference between visual and auditory tuned fraction:s
%Multi-level statistics: 
Y               = [fracFreqSes; fracOriSes];
X_mod           = [ones(nSessions,1); ones(nSessions,1)*2];
tbl             = table(Y,X_mod,[sessionData.mousename; sessionData.mousename],'VariableNames',{'Fraction','Modality','Mouse'}); %Create table for mixed model
lme             = fitlme(tbl,'Fraction~Modality+(1|Mouse)'); %construct linear mixed effects model with fixed effect of temporal window (NO random intercept for mice as cohorts are tested)
stats           = dataset2table(anova(lme,'DFMethod','Satterthwaite')); %Perform ANOVA on model and output as matrix
fprintf('Fraction auditory versus visually tuned significantly different: \n')
fprintf('F(%d,%2.0f)=%1.2f, p=%1.3f; \n',stats{2,3},stats{2,4},stats{2,2},stats{2,5})
