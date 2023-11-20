%% Oude Lohuis et al. 2024 Nat Neurosci
% Triple dissociation of auditory, visual, and motor processing in primary visual cortex
% MOL (C) 2023

% Script that analyzes V1 spiking across the session with and without muscimol

startover

%% Parameter settings:
params                      = params_histresponse_auV1;
params.AlignOn              = 'stimChange';         %On which timestamp to align trials to

params.area                 = 'V1';

params.Project              = 'CHDET'; %Which versions of the task to load data from
params.Experiment           = 'ChangeDetectionConflict'; %Which versions of the task to load data from
params.savedir              = 'E:\Documents\PhD\Figures\Project CHDET\Results - auV1\7MuscimolA1\V1responses';

params                      = MOL_getColors_CHDET(params);

%% Configurations from excel overview of recordings per project:
RootDir                 = 'E:\';
MainDataDir             = fullfile(RootDir,'Data',params.Project);
[~,~,RawLibrary]        = xlsread(fullfile(RootDir,'Data',params.Project,sprintf('%s_Sessions_Overview.xlsx',params.Project)));
RawLibrary              = RawLibrary(~cellfun(@isnumeric,RawLibrary(:,1)),:); %Trim the RawLibrary to entered values
[XlsRows,XlsColumns]    = size(RawLibrary);

ExpIdx                  = strcmp(RawLibrary(:,strcmp(RawLibrary(1,:),'MuscimolArea')),'A1');
ExpIdx = strcmp(RawLibrary(:,strcmp(RawLibrary(1,:),'MuscimolArea')),'A1') | ...
    strcmp(RawLibrary(:,strcmp(RawLibrary(1,:),'SalineArea')),'A1');
Sesselec                = RawLibrary(ExpIdx,strcmp(RawLibrary(1,:),'Rec_datetime'));

%% Get data:
[Data] = MOL_GetData('E:','CHDET',params.Experiment,{},Sesselec,{'sessionData' 'trialData_newtrials' 'spikeData'});
sessionData     = Data.sessionData;
trialData       = Data.trialData;
spikeData       = Data.spikeData;

%% Filter only V1 neurons:
spikeFields = fieldnames(spikeData);
area_idx = strcmp(spikeData.area,'V1');
for fld = 1:length(spikeFields)
    spikeData.(spikeFields{fld}) = spikeData.(spikeFields{fld})(area_idx);
end

%% Filter sessions based on containing neurons:
[sessionData,trialData,spikeData] = MOL_getTempPerSes(unique(spikeData.session_ID),sessionData,trialData,spikeData);

%% Remove sessions where V1 is affected by the muscimol
sesids = sessionData.session_ID(~ismember(sessionData.session_ID,{'20442021042332' '20452021042347'}));
[sessionData,trialData,spikeData] = MOL_getTempPerSes(sesids,sessionData,trialData,spikeData);

%% Save dataset:
save('Dataset5_1.mat','params','sessionData','trialData','spikeData')

%% Or start script from saved dataset:
load Dataset5_1.mat

%% Report dataset:
fprintf('Dataset: %d sessions, %d trials, %d neurons\n',length(sessionData.session_ID),length(trialData.session_ID),length(spikeData.session_ID));
nSessions = length(sessionData.session_ID);
nNeurons = length(spikeData.session_ID);


%% For all neurons compute histogram across session:
ratemat = NaN(nNeurons,1000);
resol = 60e6;
for iNeuron = 1:nNeurons
    
    %Get the relevant data for each session individually:
    [tempsessionData,temptrialData]     = MOL_getTempPerSes(spikeData.session_ID(iNeuron),sessionData,trialData);
    
%     edges                               = linspace(tempsessionData.t_start,tempsessionData.t_stop,101);
    edges                               = tempsessionData.t_start:resol:tempsessionData.t_stop;
    
    binnedsession                       = histcounts(spikeData.ts{iNeuron},edges) / (edges(2)-edges(1)) * 1e6;
    ratemat(iNeuron,1:length(edges)-1)                  = binnedsession;
end

%% normalize to first bin to see changes relative to the start:
ratemat = ratemat ./ repmat(ratemat(:,1),1,size(ratemat,2));

%% Make figure:
figure; hold all; set(gcf,'units','normalized','Position',[0.2 0.2 0.22 0.11],'color','w');

idx_sal         = ismember(spikeData.session_ID,sessionData.session_ID(strcmp(sessionData.SalineArea,'A1')));
idx_mus         = ismember(spikeData.session_ID,sessionData.session_ID(strcmp(sessionData.MuscimolArea,'A1')));

h = shadedErrorBar(1:size(ratemat,2),nanmean(ratemat(idx_sal,:),1),nanstd(ratemat(idx_sal,:),[],1) / sqrt(sum(idx_sal)),{'k'},0); delete(h.edge(1:2));
handles(1) = h.mainLine;
h = shadedErrorBar(1:size(ratemat,2),nanmean(ratemat(idx_mus,:),1),nanstd(ratemat(idx_mus,:),[],1) / sqrt(sum(idx_mus)),{'r'},0); delete(h.edge(1:2));
handles(2) = h.mainLine;
ylim([0 5])
ylabel('Norm. V1 rate')
xlabel('Time in session')
set(gca,'XTick',10:10:500,'XTickLabels',(10:10:500),'YTick',[1 3])
% set(gca,'XTick',[1 100],'XTickLabels',{'Start' 'End'},'YTick',[1 3])
legend(handles,{'Sal' 'Mus'},'Location','NorthWest'); legend boxoff

export_fig(fullfile(params.savedir,'Rateses_V1'),'-eps','-nocrop')

