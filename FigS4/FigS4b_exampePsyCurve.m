%% Oude Lohuis et al. 2023

%% MOL_Fig_Beh_ExampleSessions
params.Experiments         = 'BehaviorPsychophysics';
params.exAnimal           = '2.14';
params.exSession          = '2017-09-15_16-43-00';

%% General figure settings:
set(0,'defaultAxesFontSize',20)

[Data]                  = MOL_GetData('E:','DET',params.Experiments,params.exAnimal,params.exSession,{'sessionData' 'trialData'});
sessionData             = Data.sessionData;
trialData               = Data.trialData;
trialData               = MOL_RemoveLastnTrials(trialData,20); %Remove last 20 trials:

%% Save dataset:
save('DatasetS4_1.mat','params','sessionData','trialData')

%% Or start script from saved dataset:
load DatasetS4_1.mat

%% Fit and plot psychometric curve for this session
   
[visconditions,auconditions,FullRespMat,FullnTrialsMat]     = MOL_Psy_GetTrialMat(sessionData,trialData);

% Settings:
params.auprobepos       = 0.2;
params.auticks          = [0.3 0.6 1];
params.auticklabels     = ['Probe' num2cell(params.auticks*100)];
params.auxaxislabel     = 'Volume (dB)';
params.auystatslabel    = 'Auditory threshold (dB)';

params.visprobepos     = 0.005;
params.visticks        = [0.05 0.25 1];
params.vistickslabels  = ['Probe' num2cell(params.visticks*100)];
params.visxaxislabel   = 'Contrast (%)';
params.visystatslabel  = 'Visual threshold (%Contr)';

params.yticks          = [0 0.25 0.5 0.75 1];

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
[theta_est, theta_err, LLF, ctable_mod, ctable_fit, sivals, psyc_mate, ce] = mADC_model_fit_psyc_editMOL(ctable,normconditions,[]); %#ok<ASGLU>

%Store parameters:
%theta_est = 8 parameters: % 3 d' -related parameters and one 'c' parameter for each location
% dmax = theta(1:M); n = theta(M+1:2*M); s50 = theta(2*M+1:3*M); c = theta(3*M+1:end);
theta_est(5) = theta_est(5) * max(auconditions); %undo normalization of conditions
theta_est(6) = theta_est(6) * max(visconditions); %undo normalization of conditions
theta_err(5) = theta_err(5) * max(auconditions); %undo normalization of conditions
theta_err(6) = theta_err(6) * max(visconditions); %undo normalization of conditions

% Generate contingency table from fitted parameters:
[xvals_fit_au,xvals_fit_vis,ctable_fit_mat] = MOL_Gen2ADC_PsyCurve(theta_est,params);

figure; set(gcf,'color','w','units','normalized','Position', [0.1 0.1 .7 .7]);

% Construct x from position of the probe trial and conditions
xdata_au = [params.auprobepos auconditions];
% Construct x from position of the probe trial and conditions
xdata_vis = [params.visprobepos visconditions];

set(0,'Defaultlinelinewidth',5)
Colors_responses = {[1 0 0] [0 0 1] [0 0 0]};
hold all;

%Audio:
subplot(1,2,1); hold all;
%Get the hit rates for audio conditions:
%=full first dimension, vis=1 (no change), resp=1
datatoplot = squeeze(FullRespMat(:,1,1));
plot(xdata_au,datatoplot,'.','Color',Colors_responses{1},'LineWidth',5,'MarkerSize',50);
datatoplot = squeeze(FullRespMat(:,1,2));
plot(xdata_au,datatoplot,'.','Color',Colors_responses{2},'LineWidth',5,'MarkerSize',50);
%Plot single session fit, auditory responses first:
plot(xvals_fit_au,squeeze(ctable_fit_mat(1,1,:)),'k','LineWidth',3);
plot(xvals_fit_au,squeeze(ctable_fit_mat(1,2,:)),'k:','LineWidth',3);

%Visual:
subplot(1,2,2); hold all;
%Get the hit rates for visual conditions:
%=au=1 (no change), full visual dimension, resp=2
datatoplot = squeeze(FullRespMat(1,:,1));
plot(xdata_vis,datatoplot,'.','Color',Colors_responses{1},'LineWidth',5,'MarkerSize',50);
datatoplot = squeeze(FullRespMat(1,:,2));
plot(xdata_vis,datatoplot,'.','Color',Colors_responses{2},'LineWidth',5,'MarkerSize',50);
%Plot single session fit, visual responses first:
plot(xvals_fit_vis,squeeze(ctable_fit_mat(2,2,:)),'k','LineWidth',3);
plot(xvals_fit_vis,squeeze(ctable_fit_mat(2,1,:)),'k:','LineWidth',3);
title(sprintf('n=%d trials',length(trialData.trialNum)))
MOL_Psy2Sided_FigMakeup(params)
