function MOL_createPseudoPop(varargin)
%% Get input arguments:
if nargin==3
    sessionData     = varargin{1};
    trialData       = varargin{2};
    spikeData       = varargin{3};
else
    [Data] = MOL_GetData('E:','CHDET',{'ChangeDetectionConflict'},{'2003'},{},{'sessionData' 'trialData' 'spikeData'});
%     [Data] = MOL_GetData('E:','CHDET',{'ChangeDetectionConflict'},{'2012'},{'2018-08-14_14-30-15'},{'sessionData' 'trialData' 'spikeData'});
%     [Data] = MOL_GetData('E:','CHDET',{'ChangeDetectionConflict'},{'2003' '2009' '2010' '2011' '2012' '2013' '2019' '2020'},[],{'sessionData' 'trialData' 'spikeData'});
    sessionData     = Data.sessionData; %#ok<*NASGU>
    trialData       = Data.trialData;
    spikeData       = Data.spikeData;
    %% Remove last 20 trials:
    trialData = MOL_RemoveLastnTrials(trialData,20);
    spikeData = MOL_filterNeurons(sessionData,trialData,spikeData);
end

%% Make trial x neuron x time tensor by subsampling and aliging the same trial types between sessions (=Pseudopopulation)
% Be aware that the response time are not aligned
nNeurons                = length(spikeData.ts);

savedir                 = 'D:\Data\CHDET\FullDataSets\';
maxtrials               = 30; % MAXIMUM NUMBER OF TRIALS PER CONDITION
mintrials               = 10; % MININIMUM NUMBER OF TRIALS PER CONDITION

%Pseudopopulation includes visual, auditory and probe trials, grouped by
%amount of change and response. Responses to the incorrect modality are ignored. 
nConditions             = (2 * 2 * 2 + 3); %visual and auditory each have 2 levels x 2 outcomes, probe trials have 1 level, 3 outcomes. 
nTrials                 = nConditions * maxtrials;
trialTypes              = 'YYYYXXXXPPP';
visualOriChangeNorms    = [1 1 1 1 2 2 3 3 1 1 1];
audioFreqChangeNorms    = [2 2 3 3 1 1 1 1 1 1 1];
responses               = [1 3 1 3 2 3 2 3 1 2 3];

if ~isequal(length(trialTypes),length(visualOriChangeNorms),length(audioFreqChangeNorms),length(responses))
    error('unequal number of trials for selections');
end

%% Parameters for neuronal responses:
% params              = params_histresponse_pca(); % All time is in microseconds
params              = params_histresponse_DPSP();   %based on bins, without zscoring, without smoothing.
params.alignevent   = 'stimChange';                 %event timestamp to align all trials on
nTimebins           = length((params.t_pre:params.binsize:params.t_post))-1;
pseudoTensor        = NaN(nNeurons,nTrials,nTimebins); %init output matrix 

%% 
% trialTensor = true(nNeurons,nConditions);
trialTensor = NaN(nNeurons,nConditions);
lineLength = fprintf('Computing pseudopopulation: neuron      \n');

for iNeuron = 1:nNeurons
    fprintf(repmat('\b', 1, numel([num2str(iNeuron-1) num2str(nNeurons)])+2));
    fprintf('%d/%d\n',iNeuron,nNeurons);
    
%     for sesid = unique(sessionData.session_ID)'
    %% Get the relevant data for each session individually:
    [tempsessionData,temptrialData] = MOL_getTempPerSes(spikeData.session_ID(iNeuron),sessionData,trialData);
%     [temptrialData] = MOL_getTempPerSes(spikeData.session_ID(iNeuron),trialData);
    temptrialData       = MOL_CreateProbeTrials(tempsessionData,temptrialData);
    
    events_ts                           = temptrialData.(params.alignevent);          %Get events
    spikes_ts                           = spikeData.ts{iNeuron};                    %Get spikes for this neuron
    [edges,hist]                        = calc_psth(events_ts,spikes_ts,params);
    
    for iC = 1:nConditions
        selection = strcmp(temptrialData.trialType,trialTypes(iC)) ...
            & temptrialData.visualOriChangeNorm==visualOriChangeNorms(iC) ...
            & temptrialData.audioFreqChangeNorm==audioFreqChangeNorms(iC) ...
            & temptrialData.responses==responses(iC);
        
        if sum(selection)<mintrials
%             trialTensor(iNeuron,iC) = false;
        elseif sum(selection)>maxtrials
            idx = find(selection); subselec = idx(randperm(sum(selection),maxtrials));
            selection = false(size(selection));
            selection(subselec) = true;
        end
        trialTensor(iNeuron,iC) = sum(selection);
        
        idx = [1:sum(selection)] + (iC-1) * maxtrials;
        pseudoTensor(iNeuron,idx,:) = hist(selection,:);
    end
end

%% 
% removeNaNs = ~all(isnan(FullRespMat),2);
% allVar = whos;
% for iVar = 1:length(allVar)
%     if any(strfind(allVar(iVar).name,'idx_')) && ~strcmp(allVar(iVar).name,'trialidx_ses')
%         eval(sprintf('%s = %s(removeNaNs);',allVar(iVar).name,allVar(iVar).name))
%     end
% end
% FullRespMat = FullRespMat(removeNaNs,:);

params.trialTypes               = trialTypes;
params.visualOriChangeNorms     = visualOriChangeNorms;
params.audioFreqChangeNorms     = audioFreqChangeNorms;
params.responses                = responses;

params.vectrialTypes            = reshape(repmat(trialTypes,maxtrials,1),nTrials,1);
params.vecvisualOriChangeNorms  = reshape(repmat(visualOriChangeNorms,maxtrials,1),nTrials,1);
params.vecaudioFreqChangeNorms  = reshape(repmat(audioFreqChangeNorms,maxtrials,1),nTrials,1);
params.vecresponses             = reshape(repmat(responses,maxtrials,1),nTrials,1);

save(fullfile(savedir,'CHDET_PseudoPop.mat'),'edges','params','pseudoTensor','trialTensor','sessionData', 'trialData','spikeData','-v7.3');

end