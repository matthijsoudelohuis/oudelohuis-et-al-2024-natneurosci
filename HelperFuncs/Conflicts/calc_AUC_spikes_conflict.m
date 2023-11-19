function [AUC_ORI_V,AUC_ORI_C,pVal_ORI_V,pVal_ORI_C,respVis,respAud] = calc_AUC_spikes_conflict(params,sessionData,trialData,spikeData)
%% This function calculates the area under the receiver operating characteristic (auROC)
% To identify which variables are encoded in single neurons at which timepoints you can identify how well an
% external observer could discriminate variables from the firing rate at single time points. In other words,
% we used the ROC analysis (Green & Swets, 1966) to identify class-differential responses.
% To identify significant coding, we computed the area under the ROC curve for the firing rate distributions
% between two selections of trials (see below) versus a shuffled distribution.
% Matthijs oude Lohuis 2017-2023

%% Init output fields:
nNeurons                = length(spikeData.ts);

AUC_ORI_V               = NaN(nNeurons,1);
AUC_ORI_C               = NaN(nNeurons,1);

pVal_ORI_V              = NaN(nNeurons,1);
pVal_ORI_C              = NaN(nNeurons,1);

respVis                 = NaN(nNeurons,1);
respAud                 = NaN(nNeurons,1);

%% Main loop to get psth matrix:
lastsesid               = []; %memory var to keep track of whether neuron comes from different session and design matrix needs to be reconstructed

fprintf('Computing AUC for neuron        \n');

for iNeuron = 1:nNeurons %Loop over all neurons:
    fprintf(repmat('\b', 1, numel([num2str(iNeuron-1) num2str(nNeurons)])+2));
    fprintf('%d/%d\n',iNeuron,nNeurons);
    
    if ~strcmp(lastsesid,spikeData.session_ID(iNeuron)) %construct new predictor matrix if neuron comes from a new session:
        %Get the relevant data for each session individually:
        temptrialData        = MOL_getTempPerSes(spikeData.session_ID(iNeuron),trialData);
        lastsesid            = spikeData.session_ID(iNeuron); %save this session_ID
        
        %Get the right indices:
        splits                  = {};
        splits{1}       = strcmp(temptrialData.trialType,'X') & temptrialData.visualOriChangeNorm==3 & ismember(temptrialData.visualOriPostChangeNorm,[1 2]);
        splits{2}       = strcmp(temptrialData.trialType,'X') & temptrialData.visualOriChangeNorm==3 & ismember(temptrialData.visualOriPostChangeNorm,[3 4]);
        
        splits{3}      = strcmp(temptrialData.trialType,'C') & temptrialData.visualOriChangeNorm==3 & ismember(temptrialData.visualOriPostChangeNorm,[1 2]);
        splits{4}      = strcmp(temptrialData.trialType,'C') & temptrialData.visualOriChangeNorm==3 & ismember(temptrialData.visualOriPostChangeNorm,[3 4]);
    end
    
    %Compute histogram:
    events_ts               = temptrialData.(params.AlignOn);
    hist_mat                = calc_psth(events_ts,spikeData.ts{iNeuron},params);    %Construct histogram matrix:
    
    %Compute significantly responding neurons:
    bsl                     = nanmean(hist_mat(:,params.xtime>params.twin_baseline_start & params.xtime<params.twin_baseline_stop),2);
    resp                    = nanmean(hist_mat(:,params.xtime>params.twin_resp_start & params.xtime<params.twin_resp_stop),2);
    
    splitsign               = false(4,1);
    for iSplit = 1:4 %Store the mean response for each of these splits
        if sum(splits{iSplit})>=params.min_ntrial
            [~,splitsign(iSplit)]   = signrank(bsl(splits{iSplit}),resp(splits{iSplit}),'alpha',params.alpha);
        end
    end
    
    respVis(iNeuron,1)       = any(splitsign(1:2));
    respAud(iNeuron,1)        = any(splitsign(3:4));
    
    %Orientation coding: 
    feature     = [ones(sum(splits{1}),1); ones(sum(splits{2}),1)*2]; %Make correct feature vector (ones for condition A and two's for condition B)
    timeidx     = params.xtime>params.twin_resp_start & params.xtime<=params.twin_resp_stop;
    responses   = [nanmean(hist_mat(splits{1},timeidx),2); nanmean(hist_mat(splits{2},timeidx),2)]; %get all firing rates in this bin for all trials
    
    if sum(splits{1})>params.min_ntrial && sum(splits{2})>params.min_ntrial %for both conditions the number of trials is met
        
        [~,~,~,AUC_ORI_V(iNeuron,1)]          = perfcurve2(feature,responses,1); %compute AUC
        AUC_ORI_V(iNeuron,1)                  = AUC_ORI_V(iNeuron,1)*2-1;  %-1 to 1

        %Make shuffle distribution and compute permutation test threshold:
        AUC_shuffle = NaN(1,params.nshuffle); %init vector for storing AUC shuffled values
        parfor iShuf = 1:params.nshuffle
            [~,~,~,AUC_shuffle(iShuf)]    = perfcurve2(feature(randperm(length(feature))),responses,1); %#ok<PFBNS>
        end
        
        AUC_shuffle             = AUC_shuffle*2-1;
        pVal_ORI_V(iNeuron,1)     = sum(AUC_shuffle>AUC_ORI_V(iNeuron,1))/params.nshuffle;
    end
    
    %Frequency coding: 
    feature     = [ones(sum(splits{3}),1); ones(sum(splits{4}),1)*2]; %Make correct feature vector (ones for condition A and two's for condition B)
    timeidx     = params.xtime>params.twin_resp_start & params.xtime<=params.twin_resp_stop;
    responses   = [nanmean(hist_mat(splits{3},timeidx),2); nanmean(hist_mat(splits{4},timeidx),2)]; %get all firing rates in this bin for all trials
    
    if sum(splits{3})>params.min_ntrial && sum(splits{4})>params.min_ntrial %for both conditions the number of trials is met
        
        [~,~,~,AUC_ORI_C(iNeuron,1)]          = perfcurve2(feature,responses,1); %compute AUC
        AUC_ORI_C(iNeuron,1)                  = AUC_ORI_C(iNeuron,1)*2-1;  %-1 to 1
        
        %Make shuffle distribution and compute permutation test threshold:
        AUC_shuffle = NaN(1,params.nshuffle); %init vector for storing AUC shuffled values
        parfor iShuf = 1:params.nshuffle
            [~,~,~,AUC_shuffle(iShuf)]    = perfcurve2(feature(randperm(length(feature))),responses,1); %#ok<PFBNS>
        end
        
        AUC_shuffle             = AUC_shuffle*2-1;
        pVal_ORI_C(iNeuron,1)    = sum(AUC_shuffle>AUC_ORI_C(iNeuron,1))/params.nshuffle;
    end
    
end

pVal_ORI_V(pVal_ORI_V==0) = 1/params.nshuffle;
pVal_ORI_C(pVal_ORI_C==0) = 1/params.nshuffle;

end

