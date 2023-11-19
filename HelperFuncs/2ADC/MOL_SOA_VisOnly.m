%% Script that analyzes primary behavioral measures of performance across the three task versions:
% MOL (C) 2020

%% 


params.Experiments          = {'VisOnlyPsychophysics'};
% params.ExperimentLabels     = {'NE' 'UST' 'MST'};

params.nExperiments         = length(params.Experiments);

params                      = MOL_getColors_CHDET(params);

% params.savedir              = 'E:\Documents\PhD\Figures\Project CHDET\Results - 2nd Bump\1 Behavior\Supp Fig 1 - Behavior Hz Oct\';

%% General figure settings:
set(0,'defaultAxesFontSize',20)
set(0,'Defaultlinelinewidth',5)
set(0,'DefaultAxesFontName','Arial')

%% Get data:
[Data] = MOL_GetData('E:','CHDET',params.Experiments,{},[],{'sessionData' 'trialData'});
sessionData     = Data.sessionData;
trialData       = Data.trialData;

%% Remove last 20 trials:
trialData = MOL_RemoveLastnTrials(trialData,20);

%% 

binres = 50e3;

SOAbins = -0.5e6:binres:0.5e6;
SOAs = SOAbins(1:end-1)+binres/2;

nSOAs = length(SOAs);

nTrials = length(trialData.session_ID);

nConds = 5;

respmat = NaN(nSOAs-1,nConds);
offset =  NaN(nTrials,1);

for iTrial = 2:nTrials-1
    temp = trialData.stimChange(iTrial) - trialData.stimChange([iTrial-1 iTrial+1]);
    [~,i] = min(abs(temp));
    offset(iTrial) = temp(i);
end

for iSOA = 1:nSOAs
    for iCond = 1:nConds
        idx = strcmp(trialData.trialType,'X') & trialData.visualOriChangeNorm==iCond+1 & offset>=SOAbins(iSOA) & offset<=SOAbins(iSOA+1);
        respmat(iSOA,iCond) = sum(idx & trialData.vecResponse==2) / sum(idx);
%         sum(idx)
    end
end

%%

figure; set(gcf,'color','w','units','normalized','Position', [0.5 0.5 .35 .4]); hold all;

for iCond = nConds:-1:1
    plot(SOAs,respmat(:,iCond),'-','LineWidth',3,'Color',[0.2 0.2 1]*iCond/nConds)
end
ylim([0 1])
xlim([-.5e6 0.5e6])
set(gca,'XTick',[-.5e6 -.25e6 0 .25e6 .5e6],'XTickLabel',[-.5e6 -.25e6 0 .25e6 .5e6]*1e-3)
ylabel('Hit Rate')
xlabel('SOA (neg au first)')

legend(fliplr({'Imp' 'Sub' 'Thr' 'Sup' 'Max'})); legend boxoff;




