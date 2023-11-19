%% Oude Lohuis et al. 2023 Nat Neurosci
% Triple dissociation of auditory, visual, and motor processing in primary visual cortex
% MOL (C) 2023


% MOL_MuscimolAuCtx_Behavior
% This scripts plot the behavioral results of the experiments with muscimol 
% or saline infusion in bilateral auditory cortex in the multisensory change 
% detection task

startover

%% Parameters:
params = MOL_getColors_CHDET;

params.colors_audio_opto  = {[0.83 0 0.16] [1 0.47 0.37] [1 0.4 0.4]};
params.colors_visual_opto   = {[0.15 0 0.58] [0 0.67 0.94] [0.4 0.4 1]};

params.savedir              = 'E:\Documents\PhD\Figures\Project CHDET\Results - auV1\7MuscimolA1\Behavior\';

%% Get Data: %Only animals that have either saline or muscimol manipulations:
[Data] = MOL_GetData('E:','CHDET',{'ChangeDetectionConflict'},{'2022' '2026' '2030' '2031' '2044' '2045'},[],{'sessionData' 'trialData'});
sessionData     = Data.sessionData;
trialData       = Data.trialData;

%% Remove last 20 trials:
trialData = MOL_RemoveLastnTrials(trialData,20);

%% Filter sessions based on behaving (not passive):
sesselec                    = sessionData.session_ID(strcmp(sessionData.State,'Behaving'));
[sessionData,trialData]     = MOL_getTempPerSes(sesselec,sessionData,trialData);

%% Filter sessions based on manipulations in A1:
sesselec                    = sessionData.session_ID(strcmp(sessionData.MuscimolArea,'A1') | strcmp(sessionData.SalineArea,'A1'));
[sessionData,trialData]     = MOL_getTempPerSes(sesselec,sessionData,trialData);

%% Remove sessions where V1 is affected by the muscimol
sesids = sessionData.session_ID(~ismember(sessionData.session_ID,{'20442021042332' '20452021042347'}));
[sessionData,trialData] = MOL_getTempPerSes(sesids,sessionData,trialData);

%% Correct vecOctChange field in subset of old sessions:
nSessions = length(sessionData.session_ID);
for iSes = 1:nSessions
    sessionData.vecOctChange{iSes} = sessionData.vecOctChange{iSes}(sessionData.vecOctChange{iSes}~=0);
end

%% Save dataset:
save('Dataset5_2.mat','params','sessionData','trialData')

%% Or start script from saved dataset:
load Dataset5_2.mat



%% Initialize structure for saving output fit parameters:
nSessions               = length(sessionData.session_ID);
dVis                    = NaN(nSessions,2); %Init matrix for storing all visual dprime data
dAud                    = NaN(nSessions,2); %Init matrix for storing all audio dprime data
TotalResp               = NaN(3,3,3,nSessions); %Init matrix for storing all hitrate data

cVis                    = NaN(nSessions,2); %Init matrix for storing all visual dprime data
cAud                    = NaN(nSessions,2); %Init matrix for storing all audio dprime data

%% Loop over sessions:
for iSes = 1:nSessions
    fprintf('\n\n Fitting session %d/%d\n\n',iSes,nSessions);
    sesid                                                       = sessionData.session_ID(iSes);
    [tempsessionData,temptrialData]                             = MOL_getTempPerSes(sesid,sessionData,trialData);%Get the sessionData for each session individually:
    
    trialFields = fieldnames(temptrialData);
    for iF = 1:length(trialFields)
        if any(strfind(trialFields{iF},'audio')) && iscell(temptrialData.(trialFields{iF}))
            temptrialData.(trialFields{iF}) = cell2vec(temptrialData.(trialFields{iF}))';
        end
    end
    
    [visconditions,auconditions,FullRespMat,FullnTrialsMat]     = MOL_Psy_GetTrialMat(tempsessionData,temptrialData);
    %output is Au X Vis X response matrix (dimension 3 (response): layer 1 is fraction auditory, 2 visual, 3 no response)
    
    %Correction: If psychometric protocol, take intermediate values to compare with only 2 levels present
    if numel(visconditions) > 2
        FullRespMat = FullRespMat(:,[1 3 end],:,:);
        FullnTrialsMat = FullnTrialsMat(:,[1 3 end],:);
    end
    if numel(auconditions) > 2
        FullRespMat = FullRespMat([1 3 end],:,:,:);
        FullnTrialsMat = FullnTrialsMat([1 3 end],:,:);
    end
    
    %Compute d-prime for each condition:
    TotalResp(:,:,:,iSes)     = FullRespMat;
    FullRespMat_trialn          = FullRespMat .* repmat(FullnTrialsMat,1,1,3);
    
    FullRespMat(FullRespMat==0) = 1/50; %for dprime impossible to have zero or one, becomes infinite
    FullRespMat(FullRespMat==1) = 1 - 1/50;
    
    for iTrial = 1:2
        outcome         = NaN(3,3);
        outcome(1,:)    = squeeze(FullRespMat_trialn(1,1+iTrial,:)); %visual trials
        outcome(2,:)    = squeeze(FullRespMat_trialn(1+iTrial,1,:)); %audio trials
        outcome(3,:)    = squeeze(FullRespMat_trialn(1,1,:)); %probe trials
        outcome         = outcome(:,[2 1 3]); %Swap response coding, visual first.
        
%         [dVis(iSes,iTrial),dAud(iSes,iTrial),cVis(iSes,iTrial),cAud(iSes,iTrial)] = ...
%             MOL_Fit_2ADC_Full_Session(tempsessionData,temptrialData,0,outcome);
        dVis(iSes,iTrial) = norminv(FullRespMat(1,1+iTrial,2)) - norminv(FullRespMat(1,1,2));
        dAud(iSes,iTrial) = norminv(FullRespMat(1+iTrial,1,1)) - norminv(FullRespMat(1,1,1));
    end
end

fprintf('\n\n Finished fitting behavioral model.\n\n\n')

%% Get indices of sessions with saline or muscimol:
idx_sal     = strcmp(sessionData.SalineArea,'A1');
idx_musc    = strcmp(sessionData.MuscimolArea,'A1');

%% Fig 5d: show average response rate figure for saline and muscimol:
MOL_plotMuscBehavior_Rates(params,TotalResp(:,:,:,idx_sal),TotalResp(:,:,:,idx_musc),sessionData.mousename(idx_sal),sessionData.mousename(idx_musc))
subplot(1,2,1); ylim([0 0.75])
subplot(1,2,2); ylim([0 0.75])
export_fig(fullfile(params.savedir,sprintf('Mus')),'-eps','-nocrop')

%% Show average dprime figure for saline and muscimol: (not in MS)
MOL_plotMuscBehavior_Dprime(params,dVis(idx_sal,:),dAud(idx_sal,:),dVis(idx_musc,:),dAud(idx_musc,:),sessionData.mousename(idx_sal),sessionData.mousename(idx_musc))
fprintf('(n=%d,%d sessions)\n',sum(idx_sal),sum(idx_musc))






%% Dominance with muscimol in AC: (not in MS, obvious result)

uMice = unique(sessionData.mousename);
nMice = length(uMice);

SesDominanceMat_sal         = NaN(3,3,nMice); %init output var
SesDominanceMat_musc        = NaN(3,3,nMice); %init output var

for iM = 1:nMice
    %Get all sessions from this mouse together:
    mouseid         = uMice{iM};
    sesselec_sal    = unique(sessionData.session_ID(strcmp(sessionData.mousename,mouseid) & strcmp(sessionData.SalineArea,'A1')))';
    
    if ~isempty(sesselec_sal)
        %Get the sessionData for each session individually:
        [tempsessionData,temptrialData]         = MOL_getTempPerSes(sesselec_sal,sessionData,trialData);
        %Get response rates per trial type:
        [x_vis,x_au,FullRespMat,FullnTrialsMat] = MOL_Psy_GetTrialMat(tempsessionData,temptrialData); %#ok<*ASGLU>
        %Compute dominance index by taking fraction auditory / total AU + VIS
        DominanceIndexMat = (FullRespMat(:,:,1) - FullRespMat(:,:,2))./(FullRespMat(:,:,1) + FullRespMat(:,:,2));
        %Store in total dominance matrix variable:
        if FullnTrialsMat(3,3)>10 && FullnTrialsMat(2,2)>10
            SesDominanceMat_sal(:,:,iM) = DominanceIndexMat;
        end
    end
    
    sesselec_musc       = unique(sessionData.session_ID(strcmp(sessionData.mousename,mouseid) & strcmp(sessionData.MuscimolArea,'A1')))';
    
    if ~isempty(sesselec_musc)
        %Get the sessionData for each session individually:
        [tempsessionData,temptrialData]         = MOL_getTempPerSes(sesselec_musc,sessionData,trialData);
        %Get response rates per trial type:
        [x_vis,x_au,FullRespMat,FullnTrialsMat] = MOL_Psy_GetTrialMat(tempsessionData,temptrialData);
        %Compute dominance index by taking fraction auditory / total AU + VIS
        DominanceIndexMat = (FullRespMat(:,:,1) - FullRespMat(:,:,2))./(FullRespMat(:,:,1) + FullRespMat(:,:,2));
        %Store in total dominance matrix variable:
        if FullnTrialsMat(3,3)>10 && FullnTrialsMat(2,2)>10
            SesDominanceMat_musc(:,:,iM) = DominanceIndexMat;
        end
    end
end

% Figure of the average over mice:
plotDominanceMap(nanmean(SesDominanceMat_sal,3),[-0.8 0.8]);
title('Saline')

% Figure of the average over mice:
plotDominanceMap(nanmean(SesDominanceMat_musc,3),[-0.8 0.8]);
title('Muscimol')

% Compute average dominance index for matched conditions:
DomIdx_sal      = NaN(nMice,1); %Init output
DomIdx_mus      = NaN(nMice,1); %Init output

for iM = 1:nMice %loop over mice:
    %construct index to get symmetric au+vis combinations without probe trials for this mouse:
    idx             = zeros(3,3,nMice);
    idx(:,:,iM)     = eye(3);
    idx(1,1,iM)     = 0;
    %Take the mean dominance index over these values and store:
%     DomIdx_sal(iM)      = nanmean(SesDominanceMat_sal(idx==1));
%     DomIdx_mus(iM)      = nanmean(SesDominanceMat_musc(idx==1));
    DomIdx_sal(iM)      = nanmean(SesDominanceMat_sal(3,3,iM));
    DomIdx_mus(iM)      = nanmean(SesDominanceMat_musc(3,3,iM));
end

% Figure of dominance dependence on last trial:
figure; hold all; set(gcf,'units','normalized','Position',[0.7 0.4 0.1 0.35],'color','w');
temp = randn(nMice,1)/15;
scatter(temp+1,DomIdx_sal,90,'k','filled'); %scatter the individual animals
scatter(temp+2,DomIdx_mus,90,'k','filled'); %scatter the individual animals
e = errorbar(0.8,mean(DomIdx_sal),std(DomIdx_sal)/sqrt(length(DomIdx_sal)),'bo','LineWidth',5);  %plot mean + error
e = errorbar(2.2,mean(DomIdx_mus),std(DomIdx_mus)/sqrt(length(DomIdx_mus)),'bo','LineWidth',5);  %plot mean + error
plot([0.5 2.5],[0 0],'k--','LineWidth',2); %plot 0.5 dominance line for reference
p = ranksum(DomIdx_sal,DomIdx_mus);
sigstar([1 2],p)
set(gca,'XTick',[1 2],'XTickLabel',{'Saline' 'Muscimol'},'XTickLabelRotation',45) %xtick
ylabel('Dominance')
xlim([0.5 2.5])
ylim([-1 1])


%% Fig 5e: Muscimol on RT effects:
labels = {'Vthr','Vmax','Athr','Amax'};

idx_sal = ismember(trialData.session_ID,sessionData.session_ID(strcmp(sessionData.SalineArea,'A1')));
idx_mus = ismember(trialData.session_ID,sessionData.session_ID(strcmp(sessionData.MuscimolArea,'A1')));

splits{1} = idx_sal & strcmp(trialData.trialType,'X') & trialData.visualOriChangeNorm==2 & trialData.vecResponse==2;
splits{2} = idx_mus & strcmp(trialData.trialType,'X') & trialData.visualOriChangeNorm==2 & trialData.vecResponse==2;

splits{3} = idx_sal & strcmp(trialData.trialType,'X') & trialData.visualOriChangeNorm==3 & trialData.vecResponse==2;
splits{4} = idx_mus & strcmp(trialData.trialType,'X') & trialData.visualOriChangeNorm==3 & trialData.vecResponse==2;

splits{5} = idx_sal & strcmp(trialData.trialType,'Y') & trialData.audioFreqChangeNorm==2 & trialData.vecResponse==1;
splits{6} = idx_mus & strcmp(trialData.trialType,'Y') & trialData.audioFreqChangeNorm==2 & trialData.vecResponse==1;

splits{7} = idx_sal & strcmp(trialData.trialType,'Y') & trialData.audioFreqChangeNorm==3 & trialData.vecResponse==1;
splits{8} = idx_mus & strcmp(trialData.trialType,'Y') & trialData.audioFreqChangeNorm==3 & trialData.vecResponse==1;

params.nSplits = length(splits);
params.colors_splits = [params.colors_visual_opto(1:2) params.colors_visual_opto(1:2) params.colors_audio_opto(1:2) params.colors_audio_opto(1:2)];

figure; hold all; set(gcf,'units','normalized','Position',[0.2 0.4 0.2 0.23],'color','w');

for iSplit = 1:params.nSplits
    bar(iSplit,nanmean(trialData.responseLatency(splits{iSplit})),0.6,'FaceColor',params.colors_splits{iSplit});
    z = errorbar(iSplit,nanmean(trialData.responseLatency(splits{iSplit})),nanstd(trialData.responseLatency(splits{iSplit}))/sqrt(sum(splits{iSplit})),'Color','k','LineWidth',1); %#ok<*NANSTD>
end

for iSplit = 1:2:params.nSplits
    tbl             = table([trialData.responseLatency(splits{iSplit}); trialData.responseLatency(splits{iSplit+1})],[ones(sum(splits{iSplit}),1); ones(sum(splits{iSplit+1}),1)*2],[G_mou(splits{iSplit}); G_mou(splits{iSplit+1})],'VariableNames',{'RT','Treatment','Mouse'}); %Create table for mixed model
    lme             = fitlme(tbl,'RT~Treatment+(1|Mouse)'); %construct linear mixed effects model with fixed effect of treatment and random intercept for different mice
    stats           = dataset2table(anova(lme,'DFMethod','Satterthwaite')); %Perform ANOVA on model and output as matrix
    
    fprintf('Split %d vs %d (%3.1f ms to %3.1f ms, n=%d ,%d trials; \n',...
        iSplit, iSplit+1 ,nanmedian(trialData.responseLatency(splits{iSplit}))*1e-3,nanmedian(trialData.responseLatency(splits{iSplit+1}))*1e-3,sum(splits{iSplit}),sum(splits{iSplit+1}))
    fprintf('F(%d,%2.0f)=%1.2f, p=%1.3f; \n',stats{2,3},stats{2,4},stats{2,2},stats{2,5})
    p = stats{2,5};
    sigstar([0 1]+iSplit,p)
    writetable(tbl,sprintf('SourceData_Fig5e_RT_%s_Muscimol.xlsx',labels{(iSplit+1)/2}))

end

ylabel('Reaction time (ms)'); xlabel('');
xlim([0.5 params.nSplits+0.5])
set(gca,'XTick',(1:2:params.nSplits)+0.5,'XTickLabels',{'Vthr' 'Vmax' 'Athr' 'Amax'},'YTick',[0 200e3 400e3 600e3],'YTickLabels',[0 200 400 600])

filename = sprintf('Bar_RT_SalMus_trialtypes.eps');
export_fig(fullfile(params.savedir,filename),gcf);



