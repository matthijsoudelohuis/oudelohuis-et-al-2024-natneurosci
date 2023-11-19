%% Oude Lohuis et al. 2023 Nat Neurosci
% Triple dissociation of auditory, visual, and motor processing in primary visual cortex
% MOL (C) 2023

% Script analyzes the behavioral dominance in conflict trials
% MOL_ConflictDominance

params.savedir              = 'E:\Documents\PhD\Figures\Project CHDET\Results - auV1\8ConflictBehavior';

%% Get the data:
[Data] = MOL_GetData('E:','CHDET',{'BehaviorConflict'},{},[],{'sessionData' 'trialData'});
sessionData     = Data.sessionData;
trialData       = Data.trialData;

%% Filter out sessions that do not have 4 levels of visual and auditory saliency:
idx = any([cellfun(@numel,sessionData.vecFreqChange) cellfun(@numel,sessionData.vecOctChange)]==5,2) & ...
any(cellfun(@numel,sessionData.vecOriChange)==4,2);
[sessionData, trialData]        = MOL_getTempPerSes(sessionData.session_ID(idx),sessionData,trialData);

%% Save dataset:
save('Dataset6_5.mat','params','sessionData','trialData')

%% Or load dataset
load Dataset6_5.mat

%% Report dataset:
fprintf('Dataset: %d mice, %d sessions, %d trials\n',numel(unique(sessionData.mousename)),length(sessionData.session_ID),length(trialData.session_ID));

%% 

%% General settings:
params.exampleAnimals = {'2007' '2030'};


%% Loop over mice and then sessions:
sessioncounter      = 0;

uMice               = unique(sessionData.mousename)'; %Get all unique mousenames
nMice               = length(uMice); %compute length
SesDominanceMat     = NaN(5,5,nMice); %init output var
SesnTrialsMat       = NaN(5,5,nMice);
for iM = 1:nMice
    %Get all sessions from this mouse together:
    mouseid         = uMice{iM};
    sesselec        = unique(sessionData.session_ID(strcmp(sessionData.mousename,mouseid)))';
    %Get the sessionData for each session individually:
    [tempsessionData,temptrialData]         = MOL_getTempPerSes(sesselec,sessionData,trialData);
    %Get response rates per trial type:
    [x_vis,x_au,FullRespMat,FullnTrialsMat] = MOL_Psy_GetTrialMat(tempsessionData,temptrialData);
    %Compute dominance index: (A-V) / (A+V)
%     DominanceIndexMat = FullRespMat(:,:,1)./(FullRespMat(:,:,1) + FullRespMat(:,:,2));
    DominanceIndexMat = (FullRespMat(:,:,1) - FullRespMat(:,:,2)) ./(FullRespMat(:,:,1) + FullRespMat(:,:,2));
    %Store in total dominance matrix variable:
    SesDominanceMat(:,:,iM) = DominanceIndexMat;
    SesnTrialsMat(:,:,iM)   = FullnTrialsMat;
end

%% Fig 6b S10d - Make the individual dominance maps:
for i = 1:length(params.exampleAnimals)
    iM = strcmp(uMice,params.exampleAnimals(i));
    plotDominanceMap(SesDominanceMat(:,:,iM),[-1 1]);
    title(params.exampleAnimals(i))
    filename = sprintf('Conflictmap_%s.eps',params.exampleAnimals{i});
    export_fig(fullfile(params.savedir,filename),gcf)
end

%% Fig 6c : average dominance map over mice:
plotDominanceMap(nanmean(SesDominanceMat,3),[-0.8 0.8]);
title(sprintf('Average of %d mice',size(SesDominanceMat,3)))
filename = sprintf('Conflictmap_n17mice.eps');
export_fig(fullfile(params.savedir,filename),gcf)
    
%% Compute average dominance index for matched conditions:
DomIdx = NaN(nMice,1); %Init output

for iM = 1:nMice %loop over mice:
    %construct index to get symmetric au+vis combinations without probe trials for this mouse:
    idx             = zeros(5,5,nMice);
    idx(:,:,iM)     = eye(5);
    idx(1,1,iM)     = 0;
    %Take the mean dominance index over these values and store:
    DomIdx(iM)      = mean(SesDominanceMat(idx==1));
end

%% Fig 6d: dominance index overview of mice:
figure; hold all; set(gcf,'units','normalized','Position',[0.2 0.3 0.1 0.4],'color','w');
colorlist = zeros(nMice,3); %construct a list of colors for each dot
colorlist(ismember(uMice,params.exampleAnimals),:) = [[0.2 0.3 0.9]; [0.9 0.1 0.2]]; %indicate example animals from previous plot
scatter(randn(nMice,1)/15+0.9,DomIdx,90,colorlist,'filled'); %scatter the individual animals
errorbar(1.3,nanmean(DomIdx),nanstd(DomIdx)/sqrt(length(DomIdx)),'o','LineWidth',3,'Color',[0.6 0.6 0.6],'MarkerFaceColor',[0.6 0.6 0.6]);  %plot mean + error
plot([0.5 1.5],[0 0],'k--','LineWidth',2); %plot 0.5 dominance line for reference

%Figure make up:
set(gca,'XTick',1,'XTickLabel',[],'YTick',-0.8:0.4:0.8) %remove xtick
ylabel('Dominance')
p = signrank(DomIdx);
sigstar([1.3 1.35],p)
fprintf('Mean: %1.2f, Wilcoxon signed rank test, n=%d mice, p=%1.4f\n',nanmean(DomIdx),nMice,p);
xlim([0.5 1.5])
ylim([-0.8 0.8])
filename = sprintf('DI_Average_n17mice.eps');
export_fig(fullfile(params.savedir,filename),gcf)

tbl  = table(DomIdx,'VariableNames',{'Dominance Index'}); %Create table for mixed model
writetable(tbl,'SourceData_Fig6d_DominanceIndex.xlsx')

%% Control figure: correlate dominance with volume (not in MS)
sessionData.dblSoundLevel(sessionData.dblSoundLevel==0.1) = 0.6;
meanVolume = NaN(nMice,1);
for iM = 1:nMice
    meanVolume(iM) = mean(sessionData.dblSoundLevel(ismember(sessionData.mousename,uMice(iM))));
end

figure; set(gcf,'units','normalized','Position',[0.65 0.3 0.20 0.26],'color','w');
hold all; set(gcf, 'DefaultLineLineWidth', 2);
scatter(meanVolume*100,DomIdx,90,[0 0 0],'filled')
[r,p] = corr(meanVolume,DomIdx);
fprintf('Dominance does not increase with volume (r=%1.2f, p=%1.2f)\n',r,p)
xlim([40 90]); ylim([-0.3 0.9]); set(gca,'XTick',[40 65 90],'YTick',[-0.3 0 0.3 0.6 0.9])
xlabel('Sound Pressure Level (dB)')
ylabel('Dominance')
MOL_prepfigAI();
filename = sprintf('Corr_auVol_Dom.eps');
% export_fig(fullfile(params.savedir,filename),gcf)

%% Compute psychometric performance based on multi-alternative detection model:
theta_all           = NaN(nMice,8); %init output var
%Loop over mice to get dprime:
for iM = 1:nMice
    fprintf('Fitting animal  %d/%d\n',iM,nMice);
    %Get all sessions from this mouse together:
    mouseid         = uMice{iM};
    sesselec        = unique(sessionData.session_ID(strcmp(sessionData.mousename,mouseid)))';
    %Get the sessionData for each session individually:
    [tempsessionData,temptrialData]         = MOL_getTempPerSes(sesselec,sessionData,trialData);
    %Get response rates per trial type:
    [visconditions,auconditions,FullRespMat,FullnTrialsMat] = MOL_Psy_GetTrialMat(tempsessionData,temptrialData);
    
    %Rewrite fullresponse table to contingency table for model fitting (without conflict trials)
    ctable              = NaN(3,3,numel(visconditions));
    ctable(1,:,:)       = permute(FullRespMat(2:end,1,:),[3 1 2]);              %Auditory
    ctable(2,:,:)       = permute(FullRespMat(1,2:end,:),[1 3 2]);              %Visual
    ctable(3,:,:)       = repmat(permute(FullRespMat(1,1,:),[1 3 2]),1,1,numel(visconditions));    %Probe trials
    
    %Align visual and auditory intensities
    %visual and auditory conditions should be normalized such that
    %value of 1 corresponds to expected asymptotic dmax
    normconditions = [auconditions / max(auconditions); visconditions / max(visconditions)];
    
    %Fit psychometric 2ADC model:
    %     fprintf('Fitting session %d/%d, of animal %d/%d\n',ses,length(sesselec),mou,length(mouseids));
    [theta_est, theta_err, LLF, ctable_mod, ctable_fit, sivals, psyc_mate, ce] = mADC_model_fit_psyc_editMOL(ctable,normconditions,[]);
    
    %Store parameters:
    %theta_est = 8 parameters: % 3 d' -related parameters and one 'c' parameter for each location
    % dmax = theta(1:M); n = theta(M+1:2*M); s50 = theta(2*M+1:3*M); c = theta(3*M+1:end);
    theta_est(5) = theta_est(5) * max(auconditions); %undo normalization of conditions
    theta_est(6) = theta_est(6) * max(visconditions); %undo normalization of conditions
    theta_err(5) = theta_err(5) * max(auconditions); %undo normalization of conditions
    theta_err(6) = theta_err(6) * max(visconditions); %undo normalization of conditions
    
    theta_all(iM,:) = theta_est; %store in output var
end

%% Fig S10e: Control figure: correlate dominance max performance in the auditory and visual domain:
figure; set(gcf,'units','normalized','Position',[0.05 0.5 0.20 0.26],'color','w');
hold all; set(gcf, 'DefaultLineLineWidth', 2);
scatter(theta_all(:,1),DomIdx,90,[0 0 0],'filled')
[r,p] = corr(theta_all(:,1),DomIdx);
fprintf('Dominance was not correlated with auditory performance (d-prime at maximal saliency) (r=%1.2f, p=%1.2f)\n',r,p)
xlim([1 3]); ylim([-0.3 0.9]); set(gca,'XTick',[ 1 1.5 2 2.5 3],'YTick',[-0.3 0 0.3 0.6 0.9])
xlabel('Auditory D-prime')
ylabel('Dominance')
% title('Dominance - auditory performance')
MOL_prepfigAI();
filename = sprintf('Corr_Audprime_Dom.eps');
export_fig(fullfile(params.savedir,filename),gcf)

figure; set(gcf,'units','normalized','Position',[0.45 0.5 0.20 0.26],'color','w');
hold all; set(gcf, 'DefaultLineLineWidth', 2);
scatter(theta_all(:,2),DomIdx,90,[0 0 0],'filled')
[r,p] = corr(theta_all(:,2),DomIdx);
fprintf('Dominance was not correlated with visual performance (d-prime at maximal saliency)  (r=%1.2f, p=%1.2f)\n',r,p)
xlim([0.5 2.5]); ylim([-0.3 0.9]); set(gca,'XTick',[0.5 1 1.5 2 2.5 3],'YTick',[-0.3 0 0.3 0.6 0.9])
xlabel('Visual D-prime')
ylabel('Dominance')
MOL_prepfigAI();
filename = sprintf('Corr_Visdprime_Dom.eps');
export_fig(fullfile(params.savedir,filename),gcf)

%% Control figures: does dominance correlate with reaction time 
meanRT_au   = NaN(nMice,1); %compute mean reaction time for auditory trials
meanRT_vis  = NaN(nMice,1); %and for visual trials
for iM = 1:nMice
    trialidx        = ismember(trialData.session_ID,sessionData.session_ID(ismember(sessionData.mousename,uMice(iM))));
    meanRT_au(iM)   = nanmean(trialData.responseLatency(trialidx & strcmp(trialData.trialType,'Y') & trialData.correctResponse==1));
    meanRT_vis(iM)  = nanmean(trialData.responseLatency(trialidx & strcmp(trialData.trialType,'X') & trialData.correctResponse==1));
end

figure; set(gcf,'units','normalized','Position',[0.05 0.3 0.20 0.26],'color','w');
hold all; set(gcf, 'DefaultLineLineWidth', 2);
scatter(meanRT_au*1e-3,DomIdx,90,[0 0 0],'filled')
[r,p] = corr(meanRT_au,DomIdx);
fprintf('Dominance was not correlated with mean reaction time in auditory trials (r=%1.2f, p=%1.2f)\n',r,p)
xlim([200 600]); ylim([-0.3 0.9]); set(gca,'XTick',[0 200 400 600],'YTick',[-0.3 0 0.3 0.6 0.9])
xlabel('Auditory reaction time (ms)')
ylabel('Dominance')
MOL_prepfigAI();
filename = sprintf('Corr_Aulat_Dom.eps');
export_fig(fullfile(params.savedir,filename),gcf)

figure; set(gcf,'units','normalized','Position',[0.45 0.3 0.20 0.26],'color','w');
hold all; set(gcf, 'DefaultLineLineWidth', 2);
scatter(meanRT_vis*1e-3,DomIdx,90,[0 0 0],'filled')
[r,p] = corr(meanRT_vis,DomIdx);
fprintf('Dominance was not correlated with mean reaction time in visual trials (r=%1.2f, p=%1.2f)\n',r,p)
xlim([200 800]); ylim([-0.3 0.9]); set(gca,'XTick',[0 200 400 600 800],'YTick',[-0.3 0 0.3 0.6 0.9])
xlabel('Visual reaction time (ms)')
ylabel('Dominance')
MOL_prepfigAI();
filename = sprintf('Corr_Vislat_Dom.eps');
export_fig(fullfile(params.savedir,filename),gcf)

% Save dprime and rt versus dominance index as sourcedata:
tbl  = table(theta_all(:,1),theta_all(:,2),meanRT_au,meanRT_vis,DomIdx,...
    'VariableNames',{'Au Dprime','Vis Dprime','Au RT','Vis RT','Dominance Index'}); %Create table for mixed model
writetable(tbl,'SourceData_Fig10hijk_Perf_vs_DI.xlsx')

%% Fig S10e: Control figure: 
%compute for each intensity value the performance in dprime
%Then bin the dominance in auditory responses according to dprime instead of predesignated intensity values

vis_dvalues     = NaN(nMice,5); %init output var
au_dvalues      = NaN(nMice,5);

for iM = 1:nMice
    %Get all sessions from this mouse together:
    mouseid         = uMice{iM};
    sesselec        = unique(sessionData.session_ID(strcmp(sessionData.mousename,mouseid)))';
    %Get the sessionData for each session individually:
    [tempsessionData,temptrialData]         = MOL_getTempPerSes(sesselec,sessionData,trialData);
    %Get response rates per trial type:
    [visconditions,auconditions,FullRespMat,FullnTrialsMat] = MOL_Psy_GetTrialMat(tempsessionData,temptrialData);
    %add the unisensory trial types: 
    visconditions = [0 visconditions]; %#ok<AGROW> 
    auconditions = [0 auconditions]; %#ok<AGROW>
    %compute the dprime values for each condition based on the psychometric model fit: (dmax * i^n / i^n+s50^n) 
    au_dvalues(iM,:)  = theta_all(iM,1) * auconditions.^theta_all(iM,3) ./ (auconditions.^theta_all(iM,3) + theta_all(iM,5)^theta_all(iM,3));
    vis_dvalues(iM,:)  = theta_all(iM,2) * visconditions.^theta_all(iM,4) ./ (visconditions.^theta_all(iM,4) + theta_all(iM,6)^theta_all(iM,4));
end

%Bin the dprime values: 
binedges        = -0.25:0.5:2.25;
binedges(end)   = 10; %Set highest bin to inf dprime
nBins           = length(binedges)-1;

DominanceMap_dprime = NaN(nBins,nBins,nMice);%init output var
%For each mouse and for each combination of vis and au dprimes find the dominance values and store in matrix:
for iM = 1:nMice
    for iX = 1:nBins
        for iY = 1:nBins
            au_idx = au_dvalues(iM,:)>=binedges(iX) & au_dvalues(iM,:)<binedges(iX+1);
            vis_idx = vis_dvalues(iM,:)>=binedges(iY) & vis_dvalues(iM,:)<binedges(iY+1);
            if any(au_idx) && any(vis_idx) %if such a combination exists for this mouse:
                temp = SesDominanceMat(au_idx,vis_idx,iM);
                DominanceMap_dprime(iX,iY,iM) = nanmean(temp(:));
            end
        end
    end
end

%Plot average over mice:
plotDominanceMap(nanmean(DominanceMap_dprime,3),[-0.8 0.8]);
axislabels = {'0-0.5' '0.5-1' '1-1.5' '1.5-2' '2-2.5'};
set(gca,'XTickLabels',axislabels,'YTickLabels',axislabels,'FontSize',10)

filename = sprintf('Conflictmap_Dprime_n17mice.eps');
export_fig(fullfile(params.savedir,filename),gcf)

%% Compute average dominance index for matched conditions:
DomIdx = NaN(nMice,1); %Init output

for iM = 1:nMice %loop over mice:
    %construct index to get symmetric au+vis combinations without probe trials for each mouse:
    idx             = zeros(5,5,nMice);
    idx(:,:,iM)     = eye(5);
    idx(1,1,iM)     = 0;
    idx(2,2,iM)     = 0;
    %Take the mean dominance index over these values and store:
    DomIdx(iM)      = nanmean(DominanceMap_dprime(idx==1));
end

%% Fig S10f: dominance overview of mice:
figure; hold all; set(gcf,'units','normalized','Position',[0.2 0.3 0.1 0.4],'color','w');
colorlist = zeros(nMice,3); %construct a list of colors for each dot
% colorlist(ismember(uMice,params.exampleAnimals),:) = [[0.2 0.3 0.9]; [0.9 0.1 0.2]]; %indicate example animals from previous plot
scatter(randn(nMice,1)/15+0.9,DomIdx,90,colorlist,'filled'); %scatter the individual animals
e = errorbar(1.3,nanmean(DomIdx),nanstd(DomIdx)/sqrt(length(DomIdx)),'o','LineWidth',3,'Color',[0.6 0.6 0.6],'MarkerFaceColor',[0.6 0.6 0.6]);  %plot mean + error
% plot([0.5 1.5],[0.5 0.5],'k--','LineWidth',2); %plot 0.5 dominance line for reference
plot([0.5 1.5],[0 0],'k--','LineWidth',2); %plot 0.5 dominance line for reference
set(gca,'XTick',1,'XTickLabel',[],'YTick',-0.8:0.4:0.8) %remove xtick
ylabel('Dominance')
xlim([0.5 1.5])
ylim([-0.8 0.8])
p = signrank(DomIdx);
sigstar([1.3 1.35],p)
filename = sprintf('DI_Dprime_Average_n17mice.eps');
export_fig(fullfile(params.savedir,filename),gcf)
fprintf('Dominance index based on dprime different from zero, Wilcoxon signed rank test, n=%d mice, p=%1.3f\n',nMice,p);

tbl  = table(DomIdx,'VariableNames',{'Dominance Index'}); %Create table for mixed model
writetable(tbl,'SourceData_FigS10f_DominanceIndex.xlsx')

%% Fig S10g: stability of dominance between first half and second half of sessions:

SesDominanceMat_firsthalf       = NaN(5,5,nMice); %init output var
SesDominanceMat_lasthalf        = NaN(5,5,nMice); %init output var

for iM = 1:nMice
    %Get all sessions from this mouse together:
    mouseid         = uMice{iM};
    sesselec        = unique(sessionData.session_ID(strcmp(sessionData.mousename,mouseid)))';
    split = round(numel(sesselec)/2);
    sesselec_firsthalf = sesselec(1:split);
    sesselec_lasthalf = sesselec(split:end);
    
    %Get the sessionData for each session individually:
    [tempsessionData,temptrialData]         = MOL_getTempPerSes(sesselec_firsthalf,sessionData,trialData);
    %Get response rates per trial type:
    [x_vis,x_au,FullRespMat,FullnTrialsMat] = MOL_Psy_GetTrialMat(tempsessionData,temptrialData); %#ok<*ASGLU>
    %Compute dominance index by taking fraction auditory / total AU + VIS
    DominanceIndexMat = (FullRespMat(:,:,1)-FullRespMat(:,:,2))./(FullRespMat(:,:,1) + FullRespMat(:,:,2));
    %Store in total dominance matrix variable:
    SesDominanceMat_firsthalf(:,:,iM) = DominanceIndexMat;
    
    %Get the sessionData for each session individually:
    [tempsessionData,temptrialData]         = MOL_getTempPerSes(sesselec_lasthalf,sessionData,trialData);
    %Get response rates per trial type:
    [x_vis,x_au,FullRespMat,FullnTrialsMat] = MOL_Psy_GetTrialMat(tempsessionData,temptrialData);
    %Compute dominance index by taking fraction auditory / total AU + VIS
    DominanceIndexMat = (FullRespMat(:,:,1)-FullRespMat(:,:,2))./(FullRespMat(:,:,1) + FullRespMat(:,:,2));
    %Store in total dominance matrix variable:
    SesDominanceMat_lasthalf(:,:,iM) = DominanceIndexMat;
end

% Figure of the average over mice:
plotDominanceMap(nanmean(SesDominanceMat_firsthalf,3),[-0.8 0.8]);
title('First half')

% Figure of the average over mice:
plotDominanceMap(nanmean(SesDominanceMat_lasthalf,3),[-0.8 0.8]);
title('Last half')

% Compute average dominance index for matched conditions:
DomIdx_firsthalf   = NaN(nMice,1); %Init output
DomIdx_lasthalf  = NaN(nMice,1); %Init output

for iM = 1:nMice %loop over mice:
    %construct index to get symmetric au+vis combinations without probe trials for this mouse:
    idx             = zeros(5,5,nMice);
    idx(:,:,iM)     = eye(5);
    idx(1,1,iM)     = 0;
    %Take the mean dominance index over these values and store:
    DomIdx_firsthalf(iM)      = nanmean(SesDominanceMat_firsthalf(idx==1));
    DomIdx_lasthalf(iM)      = nanmean(SesDominanceMat_lasthalf(idx==1));
end

% Figure of dominance dependence on last trial:
figure; hold all; set(gcf,'units','normalized','Position',[0.2 0.3 0.15 0.6],'color','w');
temp = randn(nMice,1)/15;
scatter(temp+1.2,DomIdx_firsthalf,70,'k','filled'); %scatter the individual animals
scatter(temp+1.8,DomIdx_lasthalf,70,'k','filled'); %scatter the individual animals
for iM = 1:nMice
    plot([temp(iM)+1.2 temp(iM)+1.8],[DomIdx_firsthalf(iM) DomIdx_lasthalf(iM)],'k-','LineWidth',1); %scatter the individual animals
end
e = errorbar(0.8,mean(DomIdx_firsthalf),std(DomIdx_firsthalf)/sqrt(length(DomIdx_firsthalf)),'o','LineWidth',3,'Color',[0.6 0.6 0.6],'MarkerFaceColor',[0.6 0.6 0.6]);
e = errorbar(2.2,mean(DomIdx_lasthalf),std(DomIdx_lasthalf)/sqrt(length(DomIdx_lasthalf)),'o','LineWidth',3,'Color',[0.6 0.6 0.6],'MarkerFaceColor',[0.6 0.6 0.6]);
plot([0.5 2.5],[0 0],'k--','LineWidth',2); %plot 0.5 dominance line for reference
set(gca,'XTick',[1 2],'XTickLabel',{'First half' 'Second half'},'XTickLabelRotation',45) %xtick
ylabel('Dominance')
p = signrank(DomIdx_firsthalf,DomIdx_lasthalf);
sigstar([0.8 2.2],p)
xlim([0.5 2.5])
ylim([-1 1])
filename = sprintf('DI_lasttrial_Average_n17mice.eps');
export_fig(fullfile(params.savedir,filename),gcf)
fprintf('Dominance is stable for early and late sessions (Wilcoxon signed rank test, n=%d mice, p=%1.3f\n',nMice,p);

tbl  = table(DomIdx_firsthalf,DomIdx_lasthalf,'VariableNames',{'DI first half','DI second half'}); %Create table for mixed model
writetable(tbl,'SourceData_FigS10g_DI_sessionsplit.xlsx')

%% last trial effect (not in MS:)


%% Loop over mice and then sessions:
SesDominanceMat_lastau     = NaN(5,5,nMice); %init output var
SesDominanceMat_lastvis    = NaN(5,5,nMice); %init output var

for iM = 1:nMice
    %Get all sessions from this mouse together:
    mouseid         = uMice{iM};
    sesselec        = unique(sessionData.session_ID(strcmp(sessionData.mousename,mouseid)))';
    
    %Get the sessionData for each session individually:
    [tempsessionData,temptrialData]         = MOL_getTempPerSes(sesselec,sessionData,trialData);
    trialidx = temptrialData.correctResponse==1 & temptrialData.vecResponse==1;
    trialfields = fieldnames(temptrialData);
    for iF = 1:length(trialfields)
        temptrialData.(trialfields{iF}) = temptrialData.(trialfields{iF})([false; trialidx(1:end-1)]);
    end
    
    %Get response rates per trial type:
    [x_vis,x_au,FullRespMat,FullnTrialsMat] = MOL_Psy_GetTrialMat(tempsessionData,temptrialData);
    %Compute dominance index by taking fraction auditory / total AU + VIS
    DominanceIndexMat = FullRespMat(:,:,1)./(FullRespMat(:,:,1) + FullRespMat(:,:,2));
    %Store in total dominance matrix variable:
    SesDominanceMat_lastau(:,:,iM) = DominanceIndexMat;
    
    %Get the sessionData for each session individually:
    [tempsessionData,temptrialData]         = MOL_getTempPerSes(sesselec,sessionData,trialData);
    trialidx = temptrialData.correctResponse==1 & temptrialData.vecResponse==2;
    trialfields = fieldnames(temptrialData);
    for iF = 1:length(trialfields)
        temptrialData.(trialfields{iF}) = temptrialData.(trialfields{iF})([false; trialidx(1:end-1)]);
    end
    
    %Get response rates per trial type:
    [x_vis,x_au,FullRespMat,FullnTrialsMat] = MOL_Psy_GetTrialMat(tempsessionData,temptrialData);
    %Compute dominance index by taking fraction auditory / total AU + VIS
%     DominanceIndexMat = FullRespMat(:,:,1)./(FullRespMat(:,:,1) + FullRespMat(:,:,2));
    DominanceIndexMat = (FullRespMat(:,:,1) - FullRespMat(:,:,2))./(FullRespMat(:,:,1) + FullRespMat(:,:,2));
    %Store in total dominance matrix variable:
    SesDominanceMat_lastvis(:,:,iM) = DominanceIndexMat;
    
end

% Figure of the average over mice:
plotDominanceMap(nanmean(SesDominanceMat_lastau,3),[-1 1]);
title('Last trial: auditory correct')

% Figure of the average over mice:
plotDominanceMap(nanmean(SesDominanceMat_lastvis,3),[-1 1]);
title('Last trial: visual correct')

% Compute average dominance index for matched conditions:
DomIdx_lastau   = NaN(nMice,1); %Init output
DomIdx_lastvis  = NaN(nMice,1); %Init output

for iM = 1:nMice %loop over mice:
    %construct index to get symmetric au+vis combinations without probe trials for this mouse:
    idx             = zeros(5,5,nMice);
    idx(:,:,iM)     = eye(5);
    idx(1,1,iM)     = 0;
    %Take the mean dominance index over these values and store:
    DomIdx_lastau(iM)      = nanmean(SesDominanceMat_lastau(idx==1));
    DomIdx_lastvis(iM)      = nanmean(SesDominanceMat_lastvis(idx==1));
end

% Figure of dominance dependence on last trial:
figure; hold all; set(gcf,'units','normalized','Position',[0.2 0.3 0.15 0.6],'color','w');
temp = randn(nMice,1)/15;
scatter(temp+1.2,DomIdx_lastau,70,'k','filled'); %scatter the individual animals
scatter(temp+1.8,DomIdx_lastvis,70,'k','filled'); %scatter the individual animals
for iM = 1:nMice
    plot([temp(iM)+1.2 temp(iM)+1.8],[DomIdx_lastau(iM) DomIdx_lastvis(iM)],'k-','LineWidth',1); %scatter the individual animals
end
e = errorbar(0.8,mean(DomIdx_lastau),std(DomIdx_lastau)/sqrt(length(DomIdx_lastau)),'ro','LineWidth',3);  %plot mean + error
e = errorbar(2.2,mean(DomIdx_lastvis),std(DomIdx_lastvis)/sqrt(length(DomIdx_lastvis)),'bo','LineWidth',3);  %plot mean + error
plot([0.5 2.5],[0 0],'k--','LineWidth',2); %plot 0.5 dominance line for reference
set(gca,'XTick',[1 2],'XTickLabel',{'Last Auditory' 'Last Visual'},'XTickLabelRotation',45) %xtick
ylabel('Dominance')
p = signrank(DomIdx_lastau,DomIdx_lastvis);
sigstar([0.8 2.2],p)
xlim([0.5 2.5])
ylim([-1 1])
filename = sprintf('DI_lasttrial_Average_n17mice.eps');
export_fig(fullfile(params.savedir,filename),gcf)
fprintf('Dominance index split by previous trial (Wilcoxon signed rank test, n=%d mice, p=%1.3f\n',nMice,p);


