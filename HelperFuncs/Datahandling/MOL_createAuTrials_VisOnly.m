function [trial_out] = MOL_createAuTrials_VisOnly(sessionData,trialData)
%% This script takes as input the preprocessed Visual Only session and constructs individual trials from the auditory changes:

%% General parameters:
params.respwin              = 1e6;      %duration after stimulus change that licks are registered as correct:

%% Basis for making trialData is the stimulus trialData;
trial_out = struct();

%% Get the relevant data for each session individually:
for sesid = unique(sessionData.session_ID)'
    [tempsessionData,temptrialData]   = MOL_getTempPerSes(sesid,sessionData,trialData);
    
    [r] = cellfun(@(x) size(x,1),temptrialData.stimChangeAu);
    if any(r>1)
        temptrialData.stimChangeAu = cellfun(@transpose,temptrialData.stimChangeAu,'UniformOutput',false);
    end
    
    %Create trials to append:
    appendtrialData                     = struct();
    nAudioTrials                        = numel([temptrialData.stimChangeAu{:}]);
    appendtrialData.stimChange          = [temptrialData.stimChangeAu{:}]';
    auFields                            = fieldnames(temptrialData);
    for iF = 1:length(auFields)
        if any(strfind(auFields{iF},'audio')) && iscell(temptrialData.(auFields{iF})) && numel([temptrialData.(auFields{iF}){:}])==nAudioTrials
            appendtrialData.(auFields{iF})          = [temptrialData.(auFields{iF}){:}]';
            temptrialData.(auFields{iF})            = NaN(size(temptrialData.(auFields{iF})));
        end
    end
    appendtrialData.trialType           = repmat({'Y'},nAudioTrials,1);
    appendtrialData.session_ID          = repmat(temptrialData.session_ID(1),size(appendtrialData.stimChange));
    appendtrialData.visualOriChange     = zeros(size(appendtrialData.stimChange));
    appendtrialData.visualOriChangeNorm = ones(size(appendtrialData.stimChange));
    appendtrialData.leftCorrect         = ones(size(appendtrialData.stimChange));
    appendtrialData.rightCorrect        = ones(size(appendtrialData.stimChange));
    appendtrialData.hasvisualchange     = zeros(size(appendtrialData.stimChange));
    appendtrialData.hasaudiochange      = ones(size(appendtrialData.stimChange));
    appendtrialData.hasphotostim        = zeros(size(appendtrialData.stimChange));
    
    allLickSides                        = [temptrialData.lickSide{:}]; %all licksides as one vector
    allLickTimes                        = [temptrialData.lickTime{:}]; %all lick times as one vector
    
    appendtrialData.responseSide        = cell(size(appendtrialData.stimChange));
    appendtrialData.correctResponse     = zeros(size(appendtrialData.stimChange));
    appendtrialData.responseLatency     = NaN(size(appendtrialData.stimChange));
    %Loop over trials and register if licks happened to fall within correct
    %response window: should be random since decorrelated version
    for tr = 1:length(appendtrialData.stimChange)
        %Get licks that were in the responsewindow, then take first lick as response lick
        responselicks                   = allLickTimes(allLickTimes>appendtrialData.stimChange(tr,1) & allLickTimes<appendtrialData.stimChange(tr,1)+params.respwin);
        responsesides                   = allLickSides(allLickTimes>appendtrialData.stimChange(tr,1) & allLickTimes<appendtrialData.stimChange(tr,1)+params.respwin);
        if responselicks
            appendtrialData.noResponse(tr,1)        = 0;
            appendtrialData.responseSide{tr,1}      = responsesides(1);
            appendtrialData.responseLatency(tr,1)   = responselicks(1) - appendtrialData.stimChange(tr,1);
        else
            appendtrialData.noResponse(tr,1)   = 1;
        end
        
    end
    appendtrialData.correctResponse(strcmp(appendtrialData.responseSide,'L') & appendtrialData.leftCorrect) = 1;
    appendtrialData.correctResponse(strcmp(appendtrialData.responseSide,'R') & appendtrialData.rightCorrect) = 1;
    appendtrialData.correctResponse(appendtrialData.noResponse==1 & ~appendtrialData.leftCorrect & ~appendtrialData.rightCorrect ) = 1;
    
    if tempsessionData.VisualLeftCorrectSide
        appendtrialData.vecResponse(strcmp(appendtrialData.responseSide,'R'),1) = 1; %Response as if auditory
        appendtrialData.vecResponse(strcmp(appendtrialData.responseSide,'L'),1) = 2; %Response as if visual (left is visual correct side)
    else
        appendtrialData.vecResponse(strcmp(appendtrialData.responseSide,'R'),1) = 2; %Response as if visual
        appendtrialData.vecResponse(strcmp(appendtrialData.responseSide,'L'),1) = 1; %Response as if auditory
    end
    appendtrialData.vecResponse(appendtrialData.noResponse==1,1) = 3;
    
    %% Clear old auditory fields from temptrialData:
    auFields = {'audioFreqPreChange'
        'audioFreqPostChange'
        'audioFreqChange'
        'audioOctPreChange'
        'audioOctPostChange'
        'audioOctChange'
        'audioFreqPreChangeNorm'
        'audioFreqPostChangeNorm'
        'audioFreqChangeNorm'
        'audioOctPreChangeNorm'
        'audioOctPostChangeNorm'
        'audioOctChangeNorm'};
    
    for iF = 1:length(auFields)
        temptrialData.(auFields{iF}) = NaN(length(temptrialData.trialStart),1);
    end
    
    %% Add new trials and order & integrate:
    temptrialData                       = AppendStruct(temptrialData,appendtrialData);
    
    %sort the trials;
    trialfieldnames = fieldnames(temptrialData);
    [~,idx] = sort(temptrialData.stimChange);
    for iField = 1:length(trialfieldnames)
        temptrialData.(trialfieldnames{iField}) = temptrialData.(trialfieldnames{iField})(idx);
    end
    
    temptrialData.trialNum          = transpose(1:length(temptrialData.trialNum));
    temptrialData.trialNumInv       = transpose(-length(temptrialData.trialNum)+1:1:0);
    
    %% Add orientations and frequencies:
    switch temptrialData.trialType{1}
        case  'Y'
            temptrialData.visualOriPreChange(1)     = temptrialData.visualOriPreChange(find(strcmp(temptrialData.trialType,'X'),1));
            temptrialData.visualOriPostChange(1)    = temptrialData.visualOriPreChange(find(strcmp(temptrialData.trialType,'X'),1));
        case 'X'
            temptrialData.audioFreqPreChange(1)     = temptrialData.audioFreqPreChange(find(strcmp(temptrialData.trialType,'Y'),1));
            temptrialData.audioFreqPostChange(1)    = temptrialData.audioFreqPreChange(find(strcmp(temptrialData.trialType,'Y'),1));
            temptrialData.audioOctPreChange(1)     = temptrialData.audioOctPreChange(find(strcmp(temptrialData.trialType,'Y'),1));
            temptrialData.audioOctPostChange(1)    = temptrialData.audioOctPreChange(find(strcmp(temptrialData.trialType,'Y'),1));
    end
    for iT = 2:length(temptrialData.trialStart)
        if ~temptrialData.hasvisualchange(iT) && isnan(temptrialData.visualOriPreChange(iT))
            temptrialData.visualOriPreChange(iT) = temptrialData.visualOriPostChange(iT-1);
            temptrialData.visualOriPostChange(iT) = temptrialData.visualOriPreChange(iT);
        end
        if ~temptrialData.hasaudiochange(iT) && isnan(temptrialData.audioFreqPreChange(iT))
            temptrialData.audioFreqPreChange(iT) = temptrialData.audioFreqPostChange(iT-1);
            temptrialData.audioFreqPostChange(iT) = temptrialData.audioFreqPreChange(iT);
            temptrialData.audioOctPreChange(iT) = temptrialData.audioOctPostChange(iT-1);
            temptrialData.audioOctPostChange(iT) = temptrialData.audioOctPreChange(iT);
        end
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
    
    %% Error checks:
    if any(strcmp(tempsessionData.Experiment,{'VisOnlyTwolevels'}))
        if numel(unique(temptrialData.visualOriChangeNorm))~=3 %In this protocol there should be 2 levels of visual change + probe = 3
            error('amount of changes, orientations or frequencies is too much in a recording session')
        end
        if tempsessionData.Probe1 || tempsessionData.Probe2 %In a recording session the orientations and octaves should be one of a fixed 4
            if numel(unique(temptrialData.audioFreqChangeNorm))~=3 %In this protocol there should be 2 levels of visual change + probe = 3
                error('amount of changes, orientations or frequencies is too much in a recording session')
            end
            if ~ismember(numel(unique(temptrialData.visualOriPostChangeNorm(~isnan(temptrialData.visualOriPostChangeNorm)))),[0 4])
                error('amount of changes, orientations or frequencies is too much in a recording session')
            end
            if ~ismember(numel(unique(temptrialData.audioFreqPostChangeNorm(~isnan(temptrialData.audioFreqPostChangeNorm)))),[0 4])
                error('amount of changes, orientations or frequencies is too much in a recording session')
            end
        end
    end
    
    %Check that visual and auditory trials are interleaved:
    if mean(find(strcmp(temptrialData.trialType,'Y')))<0.4*length(temptrialData.trialNum) || mean(find(strcmp(temptrialData.trialType,'Y')))>0.6*length(temptrialData.trialNum) 
        error('Auditory trials not interleaved with rest of trials in creating new trials visual only')
    end    
    
    %% Add the new trials to the output struct:
    trial_out                           = AppendStruct(trial_out,temptrialData);
end

end