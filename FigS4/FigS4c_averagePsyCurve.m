%% Audiovisual detection task 
% Oude Lohuis et al. 2023
% Triple dissociation of auditory, visual and motor processing in primary visual cortex

%% Load the data:
[Data] = MOL_GetData('E:','DET',{'BehaviorPsychophysics' 'DetectionPsychophysics'},{},[],{'sessionData' 'trialData'});
sessionData                 = Data.sessionData;
trialData                   = Data.trialData;

params                      = MOL_getColors_CHDET();

%% Remove last 20 trials
trialData               = MOL_RemoveLastnTrials(trialData,20); 

%% Remove sessions that are passive: 
sesids                      = sessionData.session_ID(strcmp(sessionData.State,'Behaving'));
fprintf('Removed %d/%d sessions with passive stimulation\n',length(sessionData.session_ID)-length(sesids),length(sessionData.session_ID));
[sessionData, trialData]        = MOL_getTempPerSes(sesids,sessionData,trialData);

%% Save dataset:
save('DatasetS4_2.mat','params','sessionData','trialData')

%% Or load dataset:
load DatasetS4_2.mat


%% Fit each session:
nSessions               = length(sessionData.mousename);

FullParam               = NaN(nSessions,8); %Init output variable

for iSes = 1:nSessions %Fit each session
    fprintf('\nFitting session %d/%d \n\n',iSes,nSessions);
    %Get the data for this session only:
    [tempsessionData,temptrialData]                 = MOL_getTempPerSes(sessionData.session_ID(iSes),sessionData,trialData);
    %Get response rates per condition:
    [visconditions,auconditions,FullRespMat,~]      = MOL_Psy_GetTrialMat(tempsessionData,temptrialData);
    %Correct dimensions for some sessions:
    if numel(visconditions)==5 && numel(auconditions)==4
        FullRespMat = [FullRespMat; FullRespMat(end,:,:)];
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
    % order is daud, dvis, naud, nvis, s50aud, s50vis, caud, vis
    % dmax = theta(1:M); n = theta(M+1:2*M); s50 = theta(2*M+1:3*M); c = theta(3*M+1:end);
    theta_est(5) = theta_est(5) * max(auconditions); %undo normalization of conditions
    theta_est(6) = theta_est(6) * max(visconditions); %undo normalization of conditions
    FullParam(iSes,:)            = theta_est; %store parameters
end

%% Report mean number of trials:
fprintf('\nAverage number of trials: %4.2f\n',nanmean(sessionData.totalTrials))
fprintf('Minimum number of trials: %4.0f\n',min(sessionData.totalTrials))
fprintf('Maximum number of trials: %4.0f\n\n',max(sessionData.totalTrials))

%% Plot average rates for each cohort:

figure; set(gcf,'color','w','units','normalized','Position', [0.1 0.1 .7 .7]); hold all;
expanimals                 = unique(sessionData.mousename);

for iAnimal = 1:length(expanimals)
    
    [tempsessionData,~]         = MOL_getTempPerSes(sessionData.session_ID(strcmp(sessionData.mousename,expanimals{iAnimal})),sessionData,trialData);
    
    idx_ses                     = strcmp(sessionData.mousename,expanimals{iAnimal});
    meantheta_est               = median(FullParam(idx_ses,:),1);
    
    % Settings:
    params.auprobepos       = 0.2;
    params.auticks          = [0.4 0.7 1 1.3];
    params.auticklabels     = ['Probe' num2cell(params.auticks*100)];
    params.auxaxislabel     = 'Volume (dB)';
    params.auystatslabel    = 'Auditory threshold (dB)';
    
    params.visprobepos     = 0.005;
    params.visticks        = [0.05 0.25 0.7];
    params.vistickslabels  = ['Probe' num2cell(params.visticks*100)];
    params.visxaxislabel   = 'Contrast (%)';
    params.visystatslabel  = 'Visual threshold (%Contr)';
    
    params.yticks          = [0 0.25 0.5 0.75 1];
    
    % Generate contingency table from fitted parameters:
    [xvals_fit_au,xvals_fit_vis,ctable_fit_mat] = MOL_Gen2ADC_PsyCurve(meantheta_est,params);
    
    %Audio:
    subplot(1,2,1); hold all;
    plot(xvals_fit_au,squeeze(ctable_fit_mat(1,1,:)),'-','Color',[1 0 0],'LineWidth',3);
    plot(xvals_fit_au,squeeze(ctable_fit_mat(1,2,:)),':','Color',[0.2 0.2 1],'LineWidth',3);
    
    %Visual:
    subplot(1,2,2); hold all;
    plot(xvals_fit_vis,squeeze(ctable_fit_mat(2,2,:)),'-','Color',[0 0 1],'LineWidth',3);
    plot(xvals_fit_vis,squeeze(ctable_fit_mat(2,1,:)),':','Color',[1 0.2 0.2],'LineWidth',3);
end

MOL_Psy2Sided_FigMakeup(params)


%% 
colors_splits   = {'r' 'b' 'r' 'b' 'r' 'b' 'r' 'b'};
labels_splits   = {'D''Aud' 'D''Vis' 'nAud' 'nVis' 'Au Thr' 'Vis Thr' 'cAud' 'cVis'};


%% FIgures not included in the manuscript:


%% Make figures of maximal d-prime:

figure; set(gcf,'color','w','units','normalized','Position', [0.05 0.36 .17 .38]); hold all;
iP = 1;
plot(1 + randn(nSessions,1)*0.05,FullParam(:,iP),'.','Color','k','MarkerSize',30,'MarkerEdgeColor',colors_splits{iP},'LineWidth',1); hold all;
e = errorbar(1.2,mean(FullParam(:,iP)),std(FullParam(:,iP))/sqrt(nSessions),'ko','LineWidth',5);  %plot mean + error
iP = 2;
plot(2 + randn(nSessions,1)*0.05,FullParam(:,iP),'.','Color','k','MarkerSize',30,'MarkerEdgeColor',colors_splits{iP},'LineWidth',1); hold all;
e = errorbar(2.2,mean(FullParam(:,iP)),std(FullParam(:,iP))/sqrt(nSessions),'ko','LineWidth',5);  %plot mean + error
ylim([0 6]);     xlim([0.7 2.3])
ylabel('D-prime')
set(gca, 'XTick', 1:2,'XTickLabels', {'Visual' 'Audio'},'XTickLabelRotation',45)
[p,h] = signrank(FullParam(:,1),FullParam(:,2),'alpha',0.05);
fprintf('(Auditory vs visual dprime, %1.1f vs %1.1f (median), p=%1.3f, Wilcoxon signed rank test)\n',median(FullParam(:,1)),median(FullParam(:,2)),p)
sigstar([1 2],p)

%% Au threshold:
figure; set(gcf,'color','w','units','normalized','Position', [0.35 0.36 .17 .38]); hold all;
iP = 5;
plot(1 + randn(nSessions,1)*0.05,FullParam(:,iP),'.','Color','k','MarkerSize',30,'MarkerEdgeColor',colors_splits{iP},'LineWidth',1); hold all;
e = errorbar(1.2,mean(FullParam(:,iP)),std(FullParam(:,iP))/sqrt(nSessions),'ko','LineWidth',5);  %plot mean + error

ylim([0.2 1]);     xlim([0.7 1.3])
ylabel('Auditory Threshold (dB)')
set(gca, 'XTick',[])
set(gca, 'YTick',[0.2 0.6 1],'YTickLabels',[0.2 0.6 1]*100)

%% Visual threshold:
figure; set(gcf,'color','w','units','normalized','Position', [0.55 0.36 .17 .38]); hold all;
iP = 6;
plot(1 + randn(nSessions,1)*0.05,FullParam(:,iP),'.','Color','k','MarkerSize',30,'MarkerEdgeColor',colors_splits{iP},'LineWidth',1); hold all;
e = errorbar(1.2,mean(FullParam(:,iP)),std(FullParam(:,iP))/sqrt(nSessions),'ko','LineWidth',5);  %plot mean + error

ylim([0 0.2]);     xlim([0.7 1.3])
ylabel('Visual Threshold (%contrast)')
set(gca, 'XTick',[])
set(gca, 'YTick',[0 0.05 0.1 0.15 0.2],'YTickLabels',[0 0.05 0.1 0.15 0.2]*100)
