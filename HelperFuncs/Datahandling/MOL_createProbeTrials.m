function trial_out = MOL_createProbeTrials(sessionData,trialData)
%% This script takes as input the preprocessed trialdata sessions and inserts probe trials that have similar response distributions as the visual and auditory trials:

%% General parameters:
params.showdistcheck        = 0;        %Plot figure of distribution response latency of probe trials versus stimulus trials to verify matching
params.nprobes              = 120;      %Number of final probe trials per session:
params.binres               = 50e3;     %for histogram showing similar response latency distribution
params.minresplat           = 100e3;    %for histogram showing similar response latency distribution
params.maxresplat           = 0.7e6;    %maximum response latency considered for probe trials

params.postPeriod           = 3e6;      % time after a previous stimulus that is NOT elibilbe for a probe trial
params.prePeriod            = 1.5e6;    % time before the next stimulus that is NOT eligible for a probe trial

% params.preNoLick            = 0.2e6;    % time before timestamp that no lick can have occurred;
% params.postNoLickRej        = 1.5e6;    % time after timestamp that no lick can occur for rejection trial;

% params.respwin              = 1e6;      %duration after stimulus change that licks are registered as correct for VisOnly and Naive:

%% Basis for making trialData is the stimulus trialData;
trial_out = struct();

%% Get the relevant data for each session individually:
for sesid = unique(sessionData.session_ID)'
    if length(sessionData.session_ID)>1
        [tempsessionData,temptrialData]   = MOL_getTempPerSes(sesid,sessionData,trialData);
    else
        tempsessionData = sessionData; temptrialData = trialData;
    end
    
    % number of old probe trials with corresponding beh. response:
    for iResp = 1:3
        oldnprobes(iResp) = sum(strcmp(temptrialData.trialType,'P') & temptrialData.vecResponse==iResp);
    end
    
    % target
    params.nprobescond = repmat(round(params.nprobes/3),3,1); %number of target probe trials per response type
    
    %% Step 1: Go over all trials and see whether probe trials can be inserted:
    
    ITIdur = diff(temptrialData.stimChange);
    
    nTrials = length(temptrialData.session_ID);

    allSides                        = [temptrialData.lickSide{:}]; %all licksides as one vector
    allTimes                        = [temptrialData.lickTime{:}]; %all lick times as one vector
    
    candidates_aud = [];
    candidates_vis = [];
    candidates_rej = [];

    for iTrial = 1:nTrials-1 %loop over all trials
        if ITIdur(iTrial)> (params.prePeriod + params.postPeriod) %if there is a long iti then this is a candidate
            %Get all the licks in this iti:
            localLicks = allTimes(allTimes>(temptrialData.stimChange(iTrial) + params.postPeriod) & allTimes<(temptrialData.stimChange(iTrial+1) - params.prePeriod));
            localSides = allSides(allTimes>(temptrialData.stimChange(iTrial) + params.postPeriod) & allTimes<(temptrialData.stimChange(iTrial+1) - params.prePeriod));
            
            if ~isempty(localLicks) %if there are licks in this window then this is a candidate for a false alarm:
                temp = diff([temptrialData.stimChange(iTrial)+params.postPeriod localLicks]) > 1e6; %take licks that are preceded with 1 second no licking
                if any(temp) %if there are such licks:
                    idx = find(temp,1); %find the idx of this first lick with preceding window of no licks:
                    if tempsessionData.VisualLeftCorrectSide %flip based on modality-side pairing:
                        switch localSides(idx)
                            case 'L'
%                                 candidates{2}(end+1) = localLicks(idx); %#ok<*AGROW>
                                candidates_vis(end+1) = localLicks(idx); %#ok<*AGROW>
                            case 'R'
                                candidates_aud(end+1) = localLicks(idx);
                        end
                    else
                        switch localSides(idx)
                            case 'R'
                                candidates_vis(end+1) = localLicks(idx);
                            case 'L'
                                candidates_aud(end+1) = localLicks(idx);
                        end
                    end
                else %if there are no licks during this long iti then this is a correct rejection candidate:
%                     candidates{3}(end+1) = nanmean(temptrialData.stimChange(iTrial:iTrial+1));
                    candidates_rej(end+1) = nanmean(temptrialData.stimChange(iTrial:iTrial+1));
                end
            end
        end
    end
    
    %% Step 2: Take number of candidates needed to sum up to target:
%     if # of probes that we want is smaller than candidates then subsample:
    if (params.nprobescond(1) - oldnprobes(1)) < length(candidates_aud)
        candidates_aud         = candidates_aud(randperm(length(candidates_aud),max([params.nprobescond(1) - oldnprobes(1) 0])));
    end
    if (params.nprobescond(2) - oldnprobes(2)) < length(candidates_vis)
        candidates_vis         = candidates_vis(randperm(length(candidates_vis),max([params.nprobescond(2) - oldnprobes(2) 0])));
    end
    if (params.nprobescond(3) - oldnprobes(3)) < length(candidates_rej)
        candidates_rej         = candidates_rej(randperm(length(candidates_rej),max([params.nprobescond(3) - oldnprobes(3) 0])));
    end
    
    %% Step 3; Create new trials based on the selected timestamps:
    iTr = 1; %for auditory false alarms:
    if params.nprobescond(iTr)
        temp                                = temptrialData.responseLatency(temptrialData.vecResponse==iTr & strcmp(temptrialData.trialType,'Y'));
        if isempty(temp); temp = rand(30,1)*1e6; end
        temp                                = temp(~isnan(temp));
        temp                                = temp(randi(length(temp),[1 length(candidates_aud)]));
    else temp = [];
    end
    appendtrialData.responseLatency     = temp(:);
    appendtrialData.stimChange          = candidates_aud(:) - temp(:);
    appendtrialData.trialType           = repmat({'P'},length(candidates_aud),1);
    appendtrialData.vecResponse         = repmat(iTr,length(candidates_aud),1);
    appendtrialData.correctResponse     = zeros(length(candidates_aud),1);
    appendtrialData.noResponse          = zeros(length(candidates_aud),1);
    if tempsessionData.VisualLeftCorrectSide
        appendtrialData.responseSide        = repmat({'R'},length(candidates_aud),1);
    else appendtrialData.responseSide        = repmat({'L'},length(candidates_aud),1);
    end
    
    iTr = 2; %for visual false alarms:
    if params.nprobescond(iTr)
        temp                                = temptrialData.responseLatency(temptrialData.vecResponse==iTr & strcmp(temptrialData.trialType,'X'));
        if isempty(temp); temp = rand(30,1)*1e6; end
        temp                                = temp(~isnan(temp));
        temp                                = temp(randi(length(temp),[1 length(candidates_vis)]));
    else temp = [];
    end
    appendtrialData.responseLatency     = [appendtrialData.responseLatency; temp(:)];
    appendtrialData.stimChange          = [appendtrialData.stimChange;      candidates_vis(:) - temp(:)];
    appendtrialData.trialType           = [appendtrialData.trialType;       repmat({'P'},length(candidates_vis),1)];
    appendtrialData.vecResponse         = [appendtrialData.vecResponse;     repmat(iTr,length(candidates_vis),1)];
    appendtrialData.correctResponse     = [appendtrialData.correctResponse; zeros(length(candidates_vis),1)];
    appendtrialData.noResponse          = [appendtrialData.noResponse;      zeros(length(candidates_vis),1)];
    if tempsessionData.VisualLeftCorrectSide
        appendtrialData.responseSide        = [appendtrialData.responseSide;    repmat({'L'},length(candidates_vis),1)];
    else appendtrialData.responseSide        = [appendtrialData.responseSide;    repmat({'R'},length(candidates_vis),1)];
    end
    
    iTr = 3; %for correct rejections:
    appendtrialData.responseLatency     = [appendtrialData.responseLatency; NaN(length(candidates_rej),1)];
    appendtrialData.stimChange          = [appendtrialData.stimChange;      candidates_rej'];
    appendtrialData.trialType           = [appendtrialData.trialType;       repmat({'P'},length(candidates_rej),1)];
    appendtrialData.vecResponse         = [appendtrialData.vecResponse;     repmat(iTr,length(candidates_rej),1)];
    appendtrialData.correctResponse     = [appendtrialData.correctResponse; ones(length(candidates_rej),1)];
    appendtrialData.noResponse          = [appendtrialData.noResponse;      ones(length(candidates_rej),1)];
    appendtrialData.responseSide        = [appendtrialData.responseSide;    repmat({[]},length(candidates_rej),1)];
    
    %For all added probe trials:
    appendtrialData.session_ID          = repmat(temptrialData.session_ID(1),size(appendtrialData.responseLatency));
    appendtrialData.visualOriChange     = zeros(size(appendtrialData.responseLatency));
    appendtrialData.audioFreqChange     = zeros(size(appendtrialData.responseLatency));
    appendtrialData.visualOriChangeNorm = ones(size(appendtrialData.responseLatency));
    appendtrialData.audioFreqChangeNorm = ones(size(appendtrialData.responseLatency));
    appendtrialData.leftCorrect         = zeros(size(appendtrialData.responseLatency));
    appendtrialData.rightCorrect        = zeros(size(appendtrialData.responseLatency));
    appendtrialData.hasvisualchange     = zeros(size(appendtrialData.responseLatency));
    appendtrialData.hasaudiochange      = zeros(size(appendtrialData.responseLatency));
    appendtrialData.hasphotostim        = zeros(size(appendtrialData.responseLatency));
    
    %% Step 4: Add the new trials and sort:
    temptrialData                       = AppendStruct(temptrialData,appendtrialData);
    
    %sort the trials;
    [~,idx] = sort(temptrialData.stimChange);
    trialfieldnames     = fieldnames(temptrialData);
    for iField = 1:length(trialfieldnames)
        temptrialData.(trialfieldnames{iField}) = temptrialData.(trialfieldnames{iField})(idx);
    end
    
%     temptrialData.trialStart(1,1)                  = tempsessionData.t_start;
    
    %% Step 5: Add orientations and frequencies: (make trialData complete)
    if strcmp(temptrialData.trialType{1},'P') %if the first trial is a probe trial, find orientation and freq of subsequent first au and vis trials:
        temptrialData.visualOriPreChange(1)     = temptrialData.visualOriPreChange(find(strcmp(temptrialData.trialType,'X'),1));
        temptrialData.visualOriPostChange(1)    = temptrialData.visualOriPreChange(find(strcmp(temptrialData.trialType,'X'),1));
        switch tempsessionData.auChangeUnit{1}
            case 'Hz'
                temptrialData.audioFreqPreChange(1)     = temptrialData.audioFreqPreChange(find(strcmp(temptrialData.trialType,'Y'),1));
                temptrialData.audioFreqPostChange(1)    = temptrialData.audioFreqPreChange(find(strcmp(temptrialData.trialType,'Y'),1));
            case 'Oct'
                temptrialData.audioOctPreChange(1)      = temptrialData.audioOctPreChange(find(strcmp(temptrialData.trialType,'Y'),1));
                temptrialData.audioOctPostChange(1)     = temptrialData.audioOctPreChange(find(strcmp(temptrialData.trialType,'Y'),1));
        end
    end
    for iTr = 2:length(temptrialData.trialStart) %loop over trials and interpolate orientation and frequency for new probe trials:
        if strcmp(temptrialData.trialType{iTr},'P') %do this only for probe trials:
            temptrialData.visualOriPreChange(iTr,1)         = temptrialData.visualOriPostChange(iTr-1,1);
            temptrialData.visualOriPostChange(iTr,1)        = temptrialData.visualOriPreChange(iTr,1);
            switch tempsessionData.auChangeUnit{1}
                case 'Hz'
                    temptrialData.audioFreqPreChange(iTr,1)         = temptrialData.audioFreqPostChange(iTr-1,1);
                    temptrialData.audioFreqPostChange(iTr,1)        = temptrialData.audioFreqPreChange(iTr,1);
                case 'Oct'
                    temptrialData.audioOctPreChange(iTr,1)         = temptrialData.audioOctPostChange(iTr-1,1);
                    temptrialData.audioOctPostChange(iTr,1)        = temptrialData.audioOctPreChange(iTr,1);
            end
            
            temptrialData.trialStart(iTr,1)                 = temptrialData.trialEnd(iTr-1,1);
            if iTr==length(temptrialData.trialType) %if last trial add 3 seconds to trial end
                temptrialData.trialEnd(iTr,1)               = temptrialData.stimChange(iTr) + 5e6;
            elseif strcmp(temptrialData.trialType{iTr+1},'P') %if next trial also a probe trial:
                temptrialData.trialEnd(iTr,1)               = mean(temptrialData.stimChange(iTr:iTr+1));
            else %if next trial is a normal pre-existing trial then stop trial at next trial start:
                temptrialData.trialEnd(iTr,1)               = temptrialData.trialStart(iTr+1,1);
            end
        end
    end
    
    if ismember(sessionData.Experiment,{'VisOnlyPsychophysics' 'VisOnlyTwolevels'})
        for iTr = 2:length(temptrialData.trialStart)-1 %loop over trials and interpolate orientation and frequency for new probe trials:
            temptrialData.trialStart(iTr,1)                 = mean(temptrialData.stimChange(iTr-1:iTr));
            temptrialData.trialEnd(iTr,1)                   = mean(temptrialData.stimChange(iTr:iTr+1));
        end
    end
    
    if ~isfield(temptrialData,'audioOctChange')
        temptrialData.audioOctChange = zeros(size(temptrialData.responseLatency));
    end
    if ~isfield(temptrialData,'audioFreqChange')
        temptrialData.audioFreqChange = zeros(size(temptrialData.responseLatency));
    end
    temptrialData.audioOctChange(~temptrialData.hasaudiochange) = 0;
    temptrialData.audioFreqChange(~temptrialData.hasaudiochange) = 0;
    
    % Normalize changes over sessions:
    normfields = {'audioFreqChange' 'audioOctChange' 'audioVolChange' 'visualOriChange' 'visualContrChange'... %amount of changes
        'audioFreqPostChange' 'audioOctPostChange' 'audioVolPostChange' 'visualOriPostChange' 'visualContrPostChange' ... %all post change fields
        'audioFreqPreChange' 'audioOctPreChange' 'audioVolPreChange' 'visualOriPreChange' 'visualContrPreChange'}; %all prechange fields
    for iN = 1:length(normfields) %for each of the above fields, if it's present, normalize to sorted absolute unique entries
        if isfield(temptrialData,normfields{iN})
            if isnumeric(temptrialData.(normfields{iN}))
                temptrialData.([normfields{iN} 'Norm']) = NaN(size(temptrialData.trialStart));
                uniqueVals = unique(abs(temptrialData.(normfields{iN})(~isnan(temptrialData.(normfields{iN})))));
                for iCh = 1:length(uniqueVals)
                    temptrialData.([normfields{iN} 'Norm'])(abs(temptrialData.(normfields{iN}))==uniqueVals(iCh)) = iCh;
                end
            elseif iscell(temptrialData.(normfields{iN})) %visonly script (multiple uniqueVals in one 'visual trial')
                temptrialData.([normfields{iN} 'Norm']) = cell(size(temptrialData.trialStart));
                uniqueVals = unique(abs([temptrialData.(normfields{iN}){:}]));
                if numel(uniqueVals)<7
                    for iCh = 1:length(uniqueVals)
                        for iT = 1:length(temptrialData.(normfields{iN}))
                            temptrialData.([normfields{iN} 'Norm']){iT}(abs(temptrialData.(normfields{iN}){iT})==uniqueVals(iCh)) = iCh;
                        end
                    end
                end
            end
        end
    end
    
    %Align these two fields so that the code doesn't depend on it:
    if strcmp(tempsessionData.auChangeUnit,'Oct')
        temptrialData.audioFreqChangeNorm       = temptrialData.audioOctChangeNorm;
        temptrialData.audioFreqPostChangeNorm   = temptrialData.audioOctPostChangeNorm;
    elseif strcmp(tempsessionData.auChangeUnit,'Hz')
        temptrialData.audioOctChangeNorm        = temptrialData.audioFreqChangeNorm;
        temptrialData.audioOctPostChangeNorm    = temptrialData.audioFreqPostChangeNorm;
    end
    
    temptrialData.trialNum          = transpose(1:length(temptrialData.trialNum)); %Add new trial number
    temptrialData.trialNumInv       = transpose(-length(temptrialData.trialNum)+1:1:0); %and the inverse
    
    % Realign start and end of trials:
    temptrialData.trialStart(2:end) = temptrialData.trialStart(2:end) + (temptrialData.trialEnd(1:end-1) - temptrialData.trialStart(2:end));
    
    %% Reinsert removed licks:
%     if exist('removedLickTimes','var')
%         for iTr = 1:length(temptrialData.stimChange)
%             temptrialData.lickTime{iTr} = [temptrialData.lickTime{iTr} removedLickTimes(removedLickTimes>temptrialData.trialStart(iTr) & removedLickTimes<temptrialData.trialEnd(iTr))];
%             temptrialData.lickSide{iTr} = [temptrialData.lickSide{iTr} removedLickSides(removedLickTimes>temptrialData.trialStart(iTr) & removedLickTimes<temptrialData.trialEnd(iTr))];
%         end
%         
%         if numel([temptrialData.lickTime{:}])<numel(allTimes)*0.99 %threshold at 1% licks lost
%             error('Many licks got lost... ')
%         end
%     end
    
    %% Compute histograms, should match
    if params.showdistcheck
        trialtypes          = {'Y'      'X'         'P'             'P'};
        vecRespTypes        = [1        2           1               2];
        colorstrialType     = {[1 0.2 0.2] [0.2 0.2 1] [0.8 0.5 0.5] [0.5 0.5 0.8]};
        
        edges = params.minresplat:params.binres:params.maxresplat;
        resplat_dist = NaN(length(trialtypes),length(edges)-1);
        for iTr = 1:length(trialtypes)
            resplat_dist(iTr,:) = histcounts(temptrialData.responseLatency(strcmp(temptrialData.trialType,trialtypes{iTr}) & temptrialData.vecResponse==vecRespTypes(iTr)),edges,'Normalization','probability');
        end
        
        figure; hold all;
        for iTr = 1:length(trialtypes)
            plot(edges(1:end-1),resplat_dist(iTr,:)+0.005*iTr,'color',colorstrialType{iTr},'LineWidth',2)
        end
    end
    
    %% Add the new trials to the output struct:
    trial_out                           = AppendStruct(trial_out,temptrialData);
    
end

end
