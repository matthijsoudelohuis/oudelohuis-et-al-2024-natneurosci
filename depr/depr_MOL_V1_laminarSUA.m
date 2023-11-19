
%% Parameter settings:
params.Experiments          = {'ChangeDetectionConflictDecor' 'VisOnlyTwolevels' 'ChangeDetectionConflict' }; %Which versions of the task to load data from
params.ExperimentLabels     = {'NE' 'UST' 'MST'}; %Which versions of the task to load data from
params.AlignOn              = 'stimChange';         %On which timestamp to align trials to
params.nExperiments         = length(params.Experiments);
params                      = MOL_getColors_CHDET(params);

% Statistics:
params.ttestalpha             = 0.025;
params.posthoctest          = 'bonferroni';

% Plotting:
params.markersize           = 50;

params.savedir              = 'E:\Documents\PhD\Figures\Project CHDET\Results - auV1\V1 Char\';

%% Get data:
% [Data] = MOL_GetData('E:','CHDET',params.Experiments,{'2003' '2011' '2012' '2030' '2031'},[],{'sessionData' 'trialData_newtrials' 'spikeData'});
[Data] = MOL_GetData('E:','CHDET',params.Experiments,{},[],{'sessionData' 'trialData_newtrials' 'spikeData'});
sessionData         = Data.sessionData;
trialData           = Data.trialData;
spikeData           = Data.spikeData;

%% Remove last 20 trials:
trialData           = MOL_RemoveLastnTrials(trialData,20);

%% Filter out neurons based on quality:
spikeData           = MOL_filterNeurons(sessionData,trialData,spikeData);

%% Filter out multisensory sessions that have too low visual or auditory performance:
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

%% Filter out trials with photostimulation in V1:
% sesids              = sessionData.session_ID(~(sessionData.UseOpto & (strcmp(sessionData.PhotostimArea,'V1') | strcmp(sessionData.PhotostimArea,'PPC'))));
sesids              = sessionData.session_ID(sessionData.UseOpto & strcmp(sessionData.PhotostimArea,'V1'));
idx                 = trialData.hasphotostim==1 & ismember(trialData.session_ID,sesids);
trialFields         = fieldnames(trialData);
fprintf('Removed %d/%d trials  with photostim in V1\n',sum(idx),length(trialData.session_ID));
for iF = 1:length(trialFields)
    trialData.(trialFields{iF}) = trialData.(trialFields{iF})(~idx);
end

%% Filter neurons in area:
idx             = ismember(spikeData.area,{'A1' 'V1'});
spikeFields     = fieldnames(spikeData);
for iF = 1:length(spikeFields)
    spikeData.(spikeFields{iF}) = spikeData.(spikeFields{iF})(idx);
end

%% Get only sessions with neurons recorded:
[sessionData,trialData,spikeData] = MOL_getTempPerSes(unique(spikeData.session_ID),sessionData,trialData,spikeData);
fprintf('Dataset: %d sessions, %d trials, %d neurons\n',length(sessionData.session_ID),length(trialData.session_ID),length(spikeData.session_ID));



%% For each neuron compute sign resp to idx:
fprintf('Computing responses for session     \n');

params.nSplits          = 4;
params.minTrialCond     = 10;

nNeurons                = length(spikeData.ts);
lastsesid               = []; %memory var to keep track of whether neuron comes from different session and design matrix needs to be reconstructed
fprintf('Computing average Z-scored response for neuron        \n');
zmat                    = NaN(nNeurons,params.nTimebins,params.nSplits);

params.zscore           = 1;
% zresp                   = NaN(nNeurons,params.nSplits);
% params.subtr_baseline   = 1;

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
    
    splits                  = {};
    splits{1}               = ismember(temptrialData.trialType,{'X'}) & temptrialData.visualOriChangeNorm==3 & ismember(temptrialData.visualOriPostChangeNorm,[1 2]);
    splits{2}               = ismember(temptrialData.trialType,{'X'}) & temptrialData.visualOriChangeNorm==3 & ismember(temptrialData.visualOriPostChangeNorm,[3 4]);
    splits{3}               = ismember(temptrialData.trialType,{'Y'}) & temptrialData.audioFreqChangeNorm==3 & ismember(temptrialData.audioFreqPostChangeNorm,[1 2]);
    splits{4}               = ismember(temptrialData.trialType,{'Y'}) & temptrialData.audioFreqChangeNorm==3 & ismember(temptrialData.audioFreqPostChangeNorm,[3 4]);
    
    bsl                     = nanmean(hist_mat(:,params.xtime>params.twin_baseline_start & params.xtime<params.twin_baseline_stop),2);
    resp                    = nanmean(hist_mat(:,params.xtime>params.twin_resp_start & params.xtime<params.twin_resp_stop),2);
    
    splitsign               = false(4,1);
    for iSplit = 1:params.nSplits %Store the mean response for each of these splits
        if sum(splits{iSplit})>=params.minTrialCond
            zmat(iNeuron,:,iSplit)  = nanmean(hist_mat(splits{iSplit},:),1);
            
            [~,splitsign(iSplit)]   = signrank(bsl(splits{iSplit}),resp(splits{iSplit}),'alpha',params.ttestalpha);
        end
    end
    
    spikeData.sign_visresponse(iNeuron,1)       = any(splitsign(1:2));
    spikeData.sign_auresponse(iNeuron,1)        = any(splitsign(3:4));
    
    spikeData.z_visresponse(iNeuron,1)          = max([nanmean(resp(splits{1})) nanmean(resp(splits{2}))]);
    spikeData.z_auresponse(iNeuron,1)           = max([nanmean(resp(splits{3})) nanmean(resp(splits{4}))]);
    
end

zmat(zmat>20) = 20;


