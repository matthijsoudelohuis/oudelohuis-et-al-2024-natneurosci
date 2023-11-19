%% Script that analyzes primary behavioral measures of performance across the three task versions:
% MOL (C) 2024

%% Set parameters of data to load:
params.Experiments          = {'ChangeDetectionConflictDecor' 'VisOnlyPsychophysics' {'BehaviorConflict' 'BehaviorPsychophysics'}};
params.ExperimentLabels     = {'NE' 'UST' 'MST'};
params.nExperiments         = length(params.Experiments);

params                      = MOL_getColors_CHDET(params);

params.savedir              = 'E:\Documents\PhD\Figures\Project CHDET\Results - auV1\1BehaviorFig';

%% Get the data:
sessionData = struct(); trialData = struct();
for iExp = 1:length(params.Experiments)
    [Data]                  = MOL_GetData('E:','CHDET',params.Experiments{iExp},{},{},{'sessionData' 'trialData'});
    sessionData             = AppendStruct(sessionData,Data.sessionData);
    trialData               = AppendStruct(trialData,Data.trialData);
end
trialData               = MOL_RemoveLastnTrials(trialData,20); %Remove last 20 trials

%% Filter out multisensory sessions that have too low visual or auditory performance:
nSessions           = length(sessionData.session_ID);
visperf             = NaN(nSessions,1);
auperf              = NaN(nSessions,1);
for iSes = 1:nSessions
    sesidx          = strcmp(trialData.session_ID,sessionData.session_ID(iSes));
    maxsal = max(trialData.visualOriChangeNorm(sesidx));
    vistrialidx     = strcmp(trialData.trialType,'X') & trialData.visualOriChangeNorm==maxsal & trialData.hasphotostim~=1 & sesidx; 
    visperf(iSes)   = sum(trialData.correctResponse(vistrialidx)) / sum(vistrialidx); 
    autrialidx      = strcmp(trialData.trialType,'Y') & trialData.audioFreqChangeNorm==maxsal & trialData.hasphotostim~=1 & sesidx;
    auperf(iSes)    = sum(trialData.correctResponse(autrialidx)) / sum(autrialidx); 
end

sesids              = sessionData.session_ID(~(ismember(sessionData.Experiment,{'BehaviorConflict' 'BehaviorPsychophysics'}) & (visperf<0.3 | auperf<0.3)));
fprintf('Removed %d/%d sessions with low behavioral accuracy\n',length(sessionData.session_ID)-length(sesids),length(sessionData.session_ID));
[sessionData, trialData]        = MOL_getTempPerSes(sesids,sessionData,trialData);

%% Save dataset:
save('Dataset1_1.mat','params','sessionData','trialData')

%% Or load dataset:
load Dataset1_1.mat

%% Report mean number of trials:
fprintf('\nAverage number of trials: %4.2f\n',nanmean(sessionData.totalTrials))
fprintf('Minimum number of trials: %4.0f\n',min(sessionData.totalTrials))
fprintf('Maximum number of trials: %4.0f\n\n',max(sessionData.totalTrials))

%% Report reaction times:
idx_MST = ismember(trialData.session_ID,sessionData.session_ID(ismember(sessionData.Experiment,{'BehaviorConflict' 'BehaviorPsychophysics'})));
idx_UST = ismember(trialData.session_ID,sessionData.session_ID(ismember(sessionData.Experiment,{'VisOnlyPsychophysics'})));

idx_A   = strcmp(trialData.trialType,'Y') & trialData.vecResponse==1 & idx_MST;
idx_V   = strcmp(trialData.trialType,'X') & trialData.vecResponse==2 & (idx_MST | idx_UST);

RA      = prctile(trialData.responseLatency(idx_A),[50 25 75]) * 1e-3;
RV      = prctile(trialData.responseLatency(idx_V),[50 25 75]) * 1e-3;

fprintf('Auditory reaction times were fast (median: %3.0f ms; IQR: %3.0f-%3.0f ms) and faster than visual reaction times (median: %3.0f ms; IQR: %3.0f-%3.0f ms).\n',RA,RV)

%% Compare reaction times using multi-level statistics: 
nSessions           = length(sessionData.session_ID);

G_mou       = cell(nSessions,1);
uMice       = unique(sessionData.mousename);
for iMouse = 1:length(uMice)
    G_mou(strcmp(sessionData.mousename,uMice{iMouse})) = uMice(iMouse);
end

for iSes = 1:nSessions
    idx_ses = strcmp(trialData.session_ID,sessionData.session_ID(iSes));
    Y_RT_V(iSes) = nanmedian(trialData.responseLatency(idx_V & idx_ses));
    Y_RT_A(iSes) = nanmedian(trialData.responseLatency(idx_A & idx_ses));
end

% Comparing RT for visual and auditory: 
Y               = [Y_RT_V'; Y_RT_A']; %combine RTs
X_mod           = [ones(nSessions,1); ones(nSessions,1)*2];
tbl             = table(Y,X_mod,[G_mou; G_mou],'VariableNames',{'RT','Mod', 'Mouse'}); %Create table for mixed model
lme             = fitlme(tbl,'RT~Mod+(1|Mouse)'); %construct linear mixed effects model with fixed effect of modality and random intercept for different mice
stats           = dataset2table(anova(lme,'DFMethod','Satterthwaite')); %Perform ANOVA on model and output as matrix
fprintf('Visual versus auditory RT: (Linear Mixed Model)\n')
fprintf('(F(%d,%2.0f) = %1.2f, p=%1.2e, ANOVA)\n',stats{2,3},stats{2,4},stats{2,2},stats{2,5})

writetable(tbl,'SourceData_Fig1b_AV_RT.xlsx')

%% Get hit rates for each session:
nSessions               = length(sessionData.mousename);

AllRespMat             = NaN(nSessions,2,6); %for each session get hit rates for visual and auditory trials (sessions by modalities (2) times saliency levels (6))

for iSes = 1:nSessions %Fit each session
    fprintf('Computing hit rates for session %d/%d \n',iSes,nSessions);
    %Get the data for this session only:
    [tempsessionData,temptrialData]                 = MOL_getTempPerSes(sessionData.session_ID(iSes),sessionData,trialData);
    %Get response rates per condition:
    [visconditions,auconditions,FullRespMat,~]      = MOL_Psy_GetTrialMat(tempsessionData,temptrialData);
    %Correct dimensions for some sessions:
    if numel(visconditions)==5 && numel(auconditions)==4
        FullRespMat = [FullRespMat; FullRespMat(end,:,:)]; %#ok<*AGROW>
        auconditions = [2/256 8/256 32/256 64/256 128/256];
    end
    if numel(visconditions)==6 && numel(auconditions)==4
        FullRespMat = [FullRespMat; FullRespMat(end-1:end,:,:)];
        auconditions = [1/256 2/256 8/256 32/256 64/256 128/256];
    end
    if numel(visconditions)==4 && numel(auconditions)==3
        FullRespMat = [FullRespMat(1:2,:,:); FullRespMat(2:end,:,:)];
        auconditions = [2/256 8/256 32/256 128/256];
    end
    if numel(visconditions)==5 && numel(auconditions)==8
        idx = [1 3 5 7 8];
        FullRespMat = FullRespMat([1 idx+1],:,:);
        auconditions = auconditions(idx);
    end
    
    if numel(auconditions)==2
        AllRespMat(iSes,1,[1 3 6]) = FullRespMat(:,1,1);
    elseif numel(auconditions)==4
        AllRespMat(iSes,1,[1 3:6]) = FullRespMat(:,1,1);
    elseif numel(auconditions)==5
        AllRespMat(iSes,1,:) = FullRespMat(:,1,1);
    elseif numel(auconditions)==6
        AllRespMat(iSes,1,:) = FullRespMat([1:4 6:7],1,1);
    else fprintf('Unknown number of au conditions, %d\n',numel(auconditions))
    end
    
    if numel(visconditions)==2
        AllRespMat(iSes,2,[1 3 6]) = FullRespMat(1,:,2);
    elseif numel(visconditions)==4
        AllRespMat(iSes,2,[1 3:6]) = FullRespMat(1,:,2);
    elseif numel(visconditions)==5
        AllRespMat(iSes,2,:) = FullRespMat(1,:,2);
    elseif numel(visconditions)==6
        AllRespMat(iSes,2,:) = FullRespMat(1,[1:4 6:7],2);
    else fprintf('Unknown number of au conditions, %d\n',numel(visconditions))
    end

end

%% Plot average rates for each cohort:
figure; set(gcf,'color','w','units','normalized','Position', [0.15 0.5 .3 .3]); hold all;
for iExp = 1:params.nExperiments
    
    idx_exp_ses                 = ismember(sessionData.Experiment,params.Experiments{iExp});
    
    params.auprobepos = 0.001;
    params.auticks = [1/256 1/64 1/16 1/4 1/2];
    %                 params.auticklabels = {'Probe' '1/256' '1/64' '1/8' '1/2'};
    params.auticklabels = {'Catch' 'Imp' 'Sub' 'Thr' 'Sup' 'Max'};
    params.auxaxislabel  = 'Delta frequency (Oct)';
    params.auystatslabel = 'Auditory threshold (partial octave)';
    %         end
    
    params.visprobepos     = 0.5;
    params.visticks        = [5 30 90];
    params.visticks        = [1 3 7 20 90];
    params.vistickslabels  = ['Probe' num2cell(params.visticks)];
    params.vistickslabels = {'Catch' 'Imp' 'Sub' 'Thr' 'Sup' 'Max'};
    params.visxaxislabel   = 'Delta orientation (Degrees)';
    params.visystatslabel  = 'Visual threshold (Degrees)';
    
    params.yticks          = [0 0.25 0.5 0.75 1];
    
    %Audio:
    subplot(1,2,1); hold all;
    xdata_au        = [params.auprobepos params.auticks];
    meantoplot      = squeeze(nanmean(AllRespMat(idx_exp_ses,1,:),1));
    errortoplot     = squeeze(nanstd(AllRespMat(idx_exp_ses,1,:),[],1))  / sqrt(sum(idx_exp_ses));
    idx             = ~isnan(meantoplot);
    xdata_au        = xdata_au(idx); meantoplot = meantoplot(idx); errortoplot = errortoplot(idx);
    h = shadedErrorBar(xdata_au,meantoplot,errortoplot,{'-','Color',params.colors_experiments{iExp},'LineWidth',3},0);
    delete(h.edge(1)); delete(h.edge(2));
    %Visual:
    subplot(1,2,2); hold all;
    xdata_vis        = [params.visprobepos params.visticks];
    meantoplot      = squeeze(nanmean(AllRespMat(idx_exp_ses,2,:),1));
    errortoplot     = squeeze(nanstd(AllRespMat(idx_exp_ses,2,:),[],1)) / sqrt(sum(idx_exp_ses));
    idx             = ~isnan(meantoplot);
    xdata_vis       = xdata_vis(idx); meantoplot = meantoplot(idx); errortoplot = errortoplot(idx);
    h = shadedErrorBar(xdata_vis,meantoplot,errortoplot,{'-','Color',params.colors_experiments{iExp},'LineWidth',3},0);
    delete(h.edge(1)); delete(h.edge(2));
    MOL_Psy2Sided_FigMakeup(params)
end
