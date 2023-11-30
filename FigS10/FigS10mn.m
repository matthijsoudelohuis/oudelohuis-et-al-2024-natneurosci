%% Oude Lohuis et al. 2023 Nat Neurosci
% Triple dissociation of auditory, visual, and motor processing in primary visual cortex
% MOL (C) 2023

%% Parameters
params.savedir              = 'E:\Documents\PhD\Figures\Project CHDET\Results - auV1\8ConflictBehavior';

%% Get data:
[Data] = MOL_GetData('E:','CHDET',{'BehaviorSOA'},{'2026' '2027' '2030' '2031'},[],{'sessionData' 'trialData'});
sessionData     = Data.sessionData;
trialData       = Data.trialData;

%% Remove last 20 trials:
trialData = MOL_RemoveLastnTrials(trialData,20);

%% Save dataset:
save('Dataset6_4.mat','params','sessionData','trialData')

%% Or load dataset
load Dataset6_4.mat


%% General settings:
showIndFig          = 0;
showMeanFig         = 1;

set(0,'defaultAxesFontSize',20);

%% Loop over mice and then sessions:
mouseids = unique(sessionData.mousename)';

%% Calculate per mouse:
outputmat = NaN(length(mouseids),7,3); 
for mou = 1:length(mouseids)
    
    %Give title based on sessiondata information:
    mouseid         = mouseids{mou};
    sesselec        = unique(sessionData.session_ID(strcmp(sessionData.mousename,mouseid)))';
    [tempsessionData,temptrialData] = MOL_getTempPerSes(sesselec,sessionData,trialData);%Get the sessionData for each session individually:
    
    conflict_idx = strcmp(temptrialData.trialType,'C');
    soaconditions = unique(temptrialData.soa(~isnan(temptrialData.soa)));
    
    for iSOA = 1:length(soaconditions)
        for iRESP = 1:3
            outputmat(mou,iSOA,iRESP) = sum(temptrialData.vecResponse==iRESP & temptrialData.soa==soaconditions(iSOA) & conflict_idx)  / sum(temptrialData.soa==soaconditions(iSOA) & conflict_idx);
        end
    end
end

%% Plot the response rates to visual, auditory lick spout, or no lick for each of the mice:
% not in MS
figure; set(gcf,'color','w','units','normalized','Position', [0.1 0.5 0.2 0.2]);
hold all;
for mou = 1:length(mouseids)
    plot(soaconditions,outputmat(mou,:,1),'r','LineWidth',2);
    plot(soaconditions,outputmat(mou,:,2),'b','LineWidth',2);
    plot(soaconditions,outputmat(mou,:,3),'k','LineWidth',2);
end
xlim([-0.3 0.3])
set(gca,'XTick',soaconditions,'XTickLabels',soaconditions*1e3,'XTickLabelRotation',45,'YTick',[0 0.25 0.5 0.75 1],'FontSize',10)
ylabel('Response Rate','FontSize',12)

%% Loop over sessions:
sesids = unique(sessionData.session_ID)';
nSessions = length(sesids);

outputmat = NaN(length(sesids),7,3); 
for iSes = 1:nSessions
    
    %Give title based on sessiondata information:
    sess_id         = sesids{iSes};
    sesselec        = unique(sessionData.session_ID(strcmp(sessionData.session_ID,sess_id)))';
    [tempsessionData,temptrialData] = MOL_getTempPerSes(sesselec,sessionData,trialData);%Get the sessionData for each session individually:
    
    conflict_idx = strcmp(temptrialData.trialType,'C');
    soaconditions = unique(temptrialData.soa(~isnan(temptrialData.soa)));
    
    for iSOA = 1:length(soaconditions)
        for iRESP = 1:3
            outputmat(iSes,iSOA,iRESP) = sum(temptrialData.vecResponse==iRESP & temptrialData.soa==soaconditions(iSOA) & conflict_idx)  / sum(temptrialData.soa==soaconditions(iSOA) & conflict_idx);
        end
    end
end

%% Fig. S10m
figure; set(gcf,'color','w','units','normalized','Position', [0.1 0.5 0.172 0.2]);
hold all;
h = shadedErrorBar(soaconditions,nanmean(outputmat(:,:,1),1),nanstd(outputmat(:,:,1),1)/sqrt(length(sesids)),{'-','LineWidth',1,'Color','r'},0);
delete(h.edge(1)); delete(h.edge(2)); 
h = shadedErrorBar(soaconditions,nanmean(outputmat(:,:,2),1),nanstd(outputmat(:,:,1),1)/sqrt(length(sesids)),{'-','LineWidth',1,'Color','b'},0);
delete(h.edge(1)); delete(h.edge(2)); 
h = shadedErrorBar(soaconditions,nanmean(outputmat(:,:,3),1),nanstd(outputmat(:,:,1),1)/sqrt(length(sesids)),{'-','LineWidth',1,'Color','k'},0);
delete(h.edge(1)); delete(h.edge(2)); 
xlim([-0.3 0.3])
set(gca,'XTick',soaconditions,'XTickLabels',soaconditions*1e3,'XTickLabelRotation',45,'YTick',[0 0.25 0.5 0.75 1],'FontSize',10)
ylabel('Response Rate','FontSize',12)
xlabel('Stimulus onset asynchrony','FontSize',12)
MOL_prepfigAI

filename = sprintf('Conflict_SOA_Hitrates.eps');
% export_fig(fullfile(params.savedir,filename),gcf);

%% Compute Dominance index:
Dominancemat = (outputmat(:,:,1) - outputmat(:,:,2)) ./ (outputmat(:,:,1) + outputmat(:,:,2));

%% Fit the data with cumulative gaussian:
%Fit the mean data:
[mu, stddev, guessrate, lapserate, curve] = MOL_FitCumGauss_SOA(soaconditions,nanmean(Dominancemat(:,:),1));

%Make bootstrap estimate based on resampling from all sessions:
for iBoot = 1:1000
    idx                 = randi(nSessions,nSessions,1);
    bootdata = nanmean(Dominancemat(idx,:),1);
    [mu(iBoot), stddev(iBoot), guessrate(iBoot), lapserate(iBoot), curve(:,:,iBoot)] = MOL_FitCumGauss_SOA(soaconditions,bootdata);
end

%% Fig S10n:
figure; set(gcf,'color','w','units','normalized','Position', [0.1 0.5 0.172 0.2]);
hold all;

%Plot raw data:
h = errorbar(soaconditions,nanmean(Dominancemat(:,:),1),nanstd(Dominancemat(:,:),1)/sqrt(length(sesids)),'-','LineWidth',1,'Color','m');
errorbar_tick(h,0.01,'units')

%Plot bootstrapped fit and 95%CI:
meantoplot = prctile(curve(:,2,:),50,3);
errortoplot = [];
errortoplot(1,:) = abs(prctile(curve(:,2,:),2.5,3) - meantoplot);
errortoplot(2,:) = prctile(curve(:,2,:),97.5,3) - meantoplot;
h = shadedErrorBar(curve(:,1,1),meantoplot',errortoplot',{'-','LineWidth',1,'Color','k'},1);
delete(h.edge(1)); delete(h.edge(2));

%Plot estimate of crossover point:
xerr = [abs(prctile(mu,2.5) - prctile(mu,50)) prctile(mu,97.5) - prctile(mu,50)];
errorbarxy(prctile(mu,50),0.3,xerr(1),xerr(2),0,0,{'k.-', 'k', 'k'})
fprintf('Median crossover point: %2.1f ms (95%% CI: %2.1f-%2.1f ms)\n',prctile(mu,[50 2.5 97.5])*1e3)

%Figure make up:
xlim([-0.31 0.31])
ylim([-0.6 0.4])
set(gca,'XTick',soaconditions,'XTickLabels',soaconditions*1e3,'XTickLabelRotation',45,'YTick',[-0.6 -0.3 0 0.3],'FontSize',10)
ylabel('Dominance Index','FontSize',12)
xlabel('Stimulus onset asynchrony','FontSize',12)
filename = sprintf('Conflict_SOA_Dominance.eps');
MOL_prepfigAI
% export_fig(fullfile(params.savedir,filename),gcf);



