%% MOL_plot_example_conflict_neurons_raster

%% Reset all
startover

%% Parameter settings:
params                      = params_histresponse_coding; %Parameters for PSTH (All time is in microseconds)

params.Experiments          = {'ChangeDetectionConflict' 'ChangeDetectionConflictDecor'}; %Which versions of the task to load data from
% params.ExperimentLabels     = {'NE' 'UST' 'MST'}; %Labels for each task version
params.nExperiments         = length(params.Experiments);

params.AlignOn              = 'stimChange';      %On which timestamp to align as t=0

% params.clipzscore           = 25;

params.conv_sigma           = 50e3;
params.conv_win             = 'gaussian';

params.area                 = 'V1'; %Filter only V1 data

params                      = MOL_getColors_CHDET(params);

% params.trialcategories      = 'ProbeVis';
params.colors_ztrials       = [params.colors_trialtypes([2 1]) fliplr(params.colors_conflict_choices)];

params.labels_ztrials       = {'V' 'A' 'C-V' 'C-A'};
params.lines_ztrials        = {'-' '-' '-' '-'};

params.lines_experiments    = {'-' '-' '-'};

params.savedir              = 'E:\Documents\PhD\Figures\Project CHDET\Results - auV1\9ConflictNeural\ExampleNeurons';


%% Get data:
% [Data] = MOL_GetData('E:','CHDET',params.Experiments,{'2012'},[],{'sessionData' 'trialData_newtrials' 'spikeData'});
[Data] = MOL_GetData('E:','CHDET',params.Experiments,{},[],{'sessionData' 'trialData_newtrials' 'spikeData'});
sessionData     = Data.sessionData;
trialData       = Data.trialData;
spikeData       = Data.spikeData;

%% Remove last 20 trials:
% trialData = MOL_RemoveLastnTrials(trialData,20);

%% Filter neurons on area:
spikeFields     = fieldnames(spikeData);
idx             = ismember(spikeData.area,params.area);
for iF = 1:length(spikeFields)
    spikeData.(spikeFields{iF}) = spikeData.(spikeFields{iF})(idx);
end
fprintf('Filtered neurons based on area\n');

%% Filter out neurons based on quality:
spikeData       = MOL_filterNeurons(sessionData,trialData,spikeData);

%% Filter out MST sessions that have too low visual or auditory performance:
nSessions           = length(sessionData.session_ID);
visperf             = NaN(nSessions,1);
auperf              = NaN(nSessions,1);
for iSes = 1:nSessions
    sesidx          = strcmp(trialData.session_ID,sessionData.session_ID(iSes));
    vistrialidx     = strcmp(trialData.trialType,'X') & trialData.visualOriChangeNorm==3 & trialData.hasphotostim~=1 & sesidx; 
    visperf(iSes)   = sum(trialData.correctResponse(vistrialidx)) / sum(vistrialidx); 
    autrialidx      = strcmp(trialData.trialType,'Y') & trialData.audioFreqChangeNorm==3 & trialData.hasphotostim~=1 & sesidx;
    auperf(iSes)    = sum(trialData.correctResponse(autrialidx)) / sum(autrialidx); 
end

sesids              = sessionData.session_ID(~(strcmp(sessionData.Experiment,'ChangeDetectionConflict') & (visperf<0.3 | auperf<0.3)));
fprintf('Removed %d/%d sessions with low behavioral accuracy\n',length(sessionData.session_ID)-length(sesids),length(sessionData.session_ID));
[sessionData, trialData, spikeData]        = MOL_getTempPerSes(sesids,sessionData,trialData,spikeData);

%% Filter out sessions with muscimol:
% sesids              = sessionData.session_ID(~(sessionData.UseOpto & (strcmp(sessionData.PhotostimArea,'V1') | strcmp(sessionData.PhotostimArea,'PPC')))...
%     & cellfun(@isempty,sessionData.MuscimolArea));
% [sessionData, trialData, spikeData]        = MOL_getTempPerSes(sesids,sessionData,trialData,spikeData);
% fprintf('Removed sessions with activity manipulations\n');

sesids              = sessionData.session_ID(cellfun(@isempty,sessionData.MuscimolArea));
fprintf('Removed %d/%d sessions  with muscimol manipulations\n',length(sessionData.session_ID)-length(sesids),length(sessionData.session_ID));
[sessionData, trialData, spikeData]        = MOL_getTempPerSes(sesids,sessionData,trialData,spikeData);

%% 
sesids              = sessionData.session_ID(strcmp(sessionData.State,'Behaving'));
fprintf('Removed %d/%d sessions that are passive\n',length(sessionData.session_ID)-length(sesids),length(sessionData.session_ID));
[sessionData, trialData, spikeData]        = MOL_getTempPerSes(sesids,sessionData,trialData,spikeData);

%% Filter out trials with photostimulation in V1:
% sesids              = sessionData.session_ID(~(sessionData.UseOpto & (strcmp(sessionData.PhotostimArea,'V1') | strcmp(sessionData.PhotostimArea,'PPC'))));
sesids              = sessionData.session_ID(sessionData.UseOpto & strcmp(sessionData.PhotostimArea,'V1'));
idx                 = trialData.hasphotostim==1 & ismember(trialData.session_ID,sesids);
trialFields         = fieldnames(trialData);
fprintf('Removed %d/%d trials  with photostim in V1\n',sum(idx),length(trialData.session_ID));
for iF = 1:length(trialFields)
    trialData.(trialFields{iF}) = trialData.(trialFields{iF})(~idx);
end

%% Filter sessions that have only one postFreq and postOri during conflict trials:
% nSessions           = length(sessionData.session_ID);
% singleStimConfl     = false(nSessions,1);
% for iSes = 1:nSessions
%     postfreqs               = trialData.audioFreqPostChangeNorm(strcmp(trialData.trialType,'C') & strcmp(trialData.session_ID,sessionData.session_ID(iSes)));
%     postoris                = trialData.visualOriPostChangeNorm(strcmp(trialData.trialType,'C') & strcmp(trialData.session_ID,sessionData.session_ID(iSes)));
%     singleStimConfl(iSes)   = numel(unique(postfreqs))==1 && numel(unique(postoris))==1;
% %     numel(unique(postfreqs))
% %     numel(unique(postoris))
% end
% [sessionData,trialData,spikeData] = MOL_getTempPerSes(sessionData.session_ID(singleStimConfl),sessionData,trialData,spikeData);

%% Get only sessions with neurons recorded:
[sessionData,trialData,spikeData] = MOL_getTempPerSes(unique(spikeData.session_ID),sessionData,trialData,spikeData);
fprintf('Dataset: %d sessions, %d trials, %d neurons\n',length(sessionData.session_ID),length(trialData.session_ID),length(spikeData.session_ID));

%% Raster plot parameters:
params.exportfig    = 1;
params.t_pre        = -0.5e6;
params.t_post       = 1e6;
params.zscore       = 0;
params.minNtrials   = 7;

params.colors_ticks         = {[0 0 0] [238 77 71] [184 6 0]; [71 107 238] [185 0 235] [246 12 127]; [0 40 184] [68 12 246] [106 0 135]};
params.colors_ticks         = cellfun(@(x) x/256,params.colors_ticks,'UniformOutput',false);

%% Show raster plot of example MST neurons conflict trials of different saliencies in single V1 neurons
cell_IDs            = {};

%visually driven example cell:
cell_IDs{end+1}     = '20302020011621436';
cell_IDs{end+1}     = '20442021042411333';
cell_IDs{end+1}     = '20302020011621288';
cell_IDs{end+1}     = '20442021042811300';

%auditory driven example cell:
cell_IDs{end+1}     = '20302020011621092';
cell_IDs{end+1}     = '20112018081011125';

%Example of visual neural dominance: 
cell_IDs{end+1}     = '20442021042811376';
cell_IDs{end+1}     = '20302020011621372'; %check orientation preference

%Example of auditory neural dominance: 
cell_IDs{end+1}     = '20442021042811372';
cell_IDs{end+1}     = '20112018081011151';
cell_IDs{end+1}     = '20302020011621369';
cell_IDs{end+1}     = '20312020012221145';
cell_IDs{end+1}     = '20112018081011106';
cell_IDs{end+1}     = '20112018081011160';
cell_IDs{end+1}     = '20302020011621289';
cell_IDs{end+1}     = '20222019062711254';

% cell_IDs     = {'20302020011621369'};

MOL_plotRaster_Conflict(sessionData,trialData,spikeData,cell_IDs,params)
% distFig('Screen','Ext')
% distFig()

%% Show raster plot of example NE neurons conflict trials of different saliencies in single V1 neurons

% cell_IDs = spikeData.cell_ID(ismember(spikeData.session_ID,sessionData.session_ID(strcmp(sessionData.Experiment,'ChangeDetectionConflictDecor'))));
cell_IDs            = {};

%visually driven example cell:
cell_IDs{end+1}     = '10122019041011396';
cell_IDs{end+1}     = '10082019031311136';
cell_IDs{end+1}     = '10082019031211260';
cell_IDs{end+1}     = '10082019030831 218';
% cell_IDs{end+1}     = '10082019031311169';
% cell_IDs{end+1}     = '10122019041011320';

%auditory driven example cell:
cell_IDs{end+1}     = '10092019031211324';
cell_IDs{end+1}     = '10092019031211313';
cell_IDs{end+1}     = '10092019030831126';
cell_IDs{end+1}     = '10092019031311154';

%Example of visual neural dominance: 
cell_IDs{end+1}     = '10092019031311203';
cell_IDs{end+1}     = '10122019041211343';
cell_IDs{end+1}     = '10082019031211237';

%Example of auditory neural dominance: 
cell_IDs{end+1}     = '10122019041011283';
cell_IDs{end+1}     = '10092019030831134';

% cell_IDs = cell_IDs(1:50);
% cell_IDs = cell_IDs(51:101);
% cell_IDs = cell_IDs(101:163);
params.exportfig    = 1;

MOL_plotRaster_Conflict(sessionData,trialData,spikeData,cell_IDs,params)
% distFig('Screen','Ext')
% distFig()











%% 
nNeurons                = length(spikeData.ts);
lastsesid               = []; %memory var to keep track of whether neuron comes from different session and design matrix needs to be reconstructed
fprintf('Computing average Z-scored response for neuron        \n');

corr_stim_mat           = NaN(nNeurons,8); %8 conditions: 4 visual stim, 4 audio stim and how they correlate with the conflict trial
corr_choice_mat         = NaN(nNeurons,8); %8 conditions: 4 visual stim, 4 audio stim and how they correlate with the conflict trial

params.minTrialCond     = 5;

params.timeselec        = params.xtime>0 & params.xtime<0.5e6;

for iNeuron = 1:nNeurons %Loop over all neurons:
    fprintf(repmat('\b', 1, numel([num2str(iNeuron-1) num2str(nNeurons)])+2));
    fprintf('%d/%d\n',iNeuron,nNeurons);
    
    if ~strcmp(lastsesid,spikeData.session_ID(iNeuron)) %construct new predictor matrix if neuron comes from a new session:
        %Get the relevant data for each session individually:
        temptrialData        = MOL_getTempPerSes(spikeData.session_ID(iNeuron),trialData);
        lastsesid            = spikeData.session_ID(iNeuron); %save this session_ID
    end
    
    %Compute histogram:
    events_ts               = temptrialData.(params.AlignOn);
    hist_mat                = calc_psth(events_ts,spikeData.ts{iNeuron},params);    %Construct histogram matrix:
    
    % STIM:
    splits                  = {};
    splits{1}               = temptrialData.visualOriChangeNorm==2 & ismember(temptrialData.visualOriPostChangeNorm,[1 2]);
    splits{2}               = temptrialData.visualOriChangeNorm==2 & ismember(temptrialData.visualOriPostChangeNorm,[3 4]);
    splits{3}               = temptrialData.visualOriChangeNorm==3 & ismember(temptrialData.visualOriPostChangeNorm,[1 2]);
    splits{4}               = temptrialData.visualOriChangeNorm==3 & ismember(temptrialData.visualOriPostChangeNorm,[3 4]);
    
    nSplits                 = length(splits);
    for iSpl = 1:nSplits
        idx_X = splits{iSpl} & strcmp(temptrialData.trialType,'X');
        idx_C = splits{iSpl} & strcmp(temptrialData.trialType,'C');
        
        if sum(idx_X)>=params.minTrialCond && sum(idx_C)>=params.minTrialCond
            corr_stim_mat(iNeuron,iSpl) = corr(nanmean(hist_mat(idx_X,params.timeselec),1)',nanmean(hist_mat(idx_C,params.timeselec),1)');
        end
    end
    
    splits                  = {};
    splits{1}               = temptrialData.audioFreqChangeNorm==2 & ismember(temptrialData.audioFreqPostChangeNorm,[1 2]);
    splits{2}               = temptrialData.audioFreqChangeNorm==2 & ismember(temptrialData.audioFreqPostChangeNorm,[3 4]);
    splits{3}               = temptrialData.audioFreqChangeNorm==3 & ismember(temptrialData.audioFreqPostChangeNorm,[1 2]);
    splits{4}               = temptrialData.audioFreqChangeNorm==3 & ismember(temptrialData.audioFreqPostChangeNorm,[3 4]);
    
    nSplits                 = length(splits);
    for iSpl = 1:nSplits
        idx_Y = splits{iSpl} & strcmp(temptrialData.trialType,'Y');
        idx_C = splits{iSpl} & strcmp(temptrialData.trialType,'C');
        
        if sum(idx_Y)>=params.minTrialCond && sum(idx_C)>=params.minTrialCond
            corr_stim_mat(iNeuron,iSpl+4) = corr(nanmean(hist_mat(idx_Y,params.timeselec),1)',nanmean(hist_mat(idx_C,params.timeselec),1)');
        end
    end
    
    % CHOICE:
    splits                  = {};
    splits{1}               = temptrialData.visualOriChangeNorm==2 & ismember(temptrialData.visualOriPostChangeNorm,[1 2]);
    splits{2}               = temptrialData.visualOriChangeNorm==2 & ismember(temptrialData.visualOriPostChangeNorm,[3 4]);
    splits{3}               = temptrialData.visualOriChangeNorm==3 & ismember(temptrialData.visualOriPostChangeNorm,[1 2]);
    splits{4}               = temptrialData.visualOriChangeNorm==3 & ismember(temptrialData.visualOriPostChangeNorm,[3 4]);
    
    nSplits                 = length(splits);
    for iSpl = 1:nSplits
        idx_X = splits{iSpl} & strcmp(temptrialData.trialType,'X') & temptrialData.vecResponse==2;
        idx_C = splits{iSpl} & strcmp(temptrialData.trialType,'C') & temptrialData.vecResponse==2;
        
        if sum(idx_X)>=params.minTrialCond && sum(idx_C)>=params.minTrialCond
            corr_choice_mat(iNeuron,iSpl) = corr(nanmean(hist_mat(idx_X,params.timeselec),1)',nanmean(hist_mat(idx_C,params.timeselec),1)');
        end
    end
    
    splits                  = {};
    splits{1}               = temptrialData.audioFreqChangeNorm==2 & ismember(temptrialData.audioFreqPostChangeNorm,[1 2]);
    splits{2}               = temptrialData.audioFreqChangeNorm==2 & ismember(temptrialData.audioFreqPostChangeNorm,[3 4]);
    splits{3}               = temptrialData.audioFreqChangeNorm==3 & ismember(temptrialData.audioFreqPostChangeNorm,[1 2]);
    splits{4}               = temptrialData.audioFreqChangeNorm==3 & ismember(temptrialData.audioFreqPostChangeNorm,[3 4]);
    
    nSplits                 = length(splits);
    for iSpl = 1:nSplits
        idx_Y = splits{iSpl} & strcmp(temptrialData.trialType,'Y') & temptrialData.vecResponse==1;
        idx_C = splits{iSpl} & strcmp(temptrialData.trialType,'C') & temptrialData.vecResponse==1;
        
        if sum(idx_Y)>=params.minTrialCond && sum(idx_C)>=params.minTrialCond
            corr_choice_mat(iNeuron,iSpl+4) = corr(nanmean(hist_mat(idx_Y,params.timeselec),1)',nanmean(hist_mat(idx_C,params.timeselec),1)');
        end
    end
end

% unique(temptrialData.visualOriPostChangeNorm(ismember(temptrialData.trialType,{'C'})))
% unique(temptrialData.audioFreqPostChangeNorm(ismember(temptrialData.trialType,{'C'})))

%%

figure; hold all; set(gcf,'units','normalized','Position',[0.12 0.5 0.25 0.35],'color','w');

meanvis = nanmean(corr_stim_mat(:,1:2),2);
meanaud = nanmean(corr_stim_mat(:,5:6),2);
scatter(meanvis,meanaud,50,'g.');

meanvis = nanmean(corr_stim_mat(:,3:4),2);
meanaud = nanmean(corr_stim_mat(:,7:8),2);
scatter(meanvis,meanaud,50,'k.');

nanidx = ~isnan(meanvis) & ~isnan(meanaud);
[r,p] = corr(meanvis(nanidx),meanaud(nanidx));

%%
figure; hold all; set(gcf,'units','normalized','Position',[0.12 0.5 0.25 0.35],'color','w');

meanvis = nanmean(corr_choice_mat(:,1:2),2);
meanaud = nanmean(corr_choice_mat(:,5:6),2);
scatter(meanvis,meanaud,50,'g.');

meanvis = nanmean(corr_choice_mat(:,3:4),2);
meanaud = nanmean(corr_choice_mat(:,7:8),2);
scatter(meanvis,meanaud,50,'k.');

nanidx = ~isnan(meanvis) & ~isnan(meanaud);
[r,p] = corr(meanvis(nanidx),meanaud(nanidx));


