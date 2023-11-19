function [AUC_ORI_V,AUC_ORI_C,AUC_ORI_C_MF,AUC_ORI_C_MS] = MOL_auROC_auV1_conflict_glm(params,sessionData,trialData,spikeData,output)
% function [AUC_ORI_V,AUC_ORI_C,pVal_ORI_V,pVal_ORI_C] = MOL_auROC_auV1_conflict_glm(params,sessionData,trialData,spikeData,output)
%% This function calculates the area under the receiver operating characteristic (auROC)
% To identify which variables are encoded in single neurons at which timepoints you can identify how well an
% external observer could discriminate variables from the firing rate at single time points. In other words,
% we used the ROC analysis (Green & Swets, 1966) to identify class-differential responses.
% To identify significant coding, we computed the area under the ROC curve for the firing rate distributions
% between two selections of trials (see below) versus a shuffled distribution.
% Matthijs oude Lohuis 2017-2021

%% Init output fields:
nNeurons                = length(spikeData.ts);

AUC_ORI_V               = NaN(nNeurons,1);
AUC_ORI_C               = NaN(nNeurons,1);

AUC_ORI_C_MF            = NaN(nNeurons,1); %full model
AUC_ORI_C_MS            = NaN(nNeurons,1); %shuffled model

% pVal_ORI_V              = NaN(nNeurons,1);
% pVal_ORI_C              = NaN(nNeurons,1);

%% Main loop to get psth matrix:
lastsesid               = []; %memory var to keep track of whether neuron comes from different session and design matrix needs to be reconstructed

fprintf('Computing AUC for neuron        \n');

%construct session-shuffled redistribution:

params.lambdastring = 'lambda_min';

output.modelFits_sh = output.modelFits;
weightmat = NaN(nNeurons,params.nPredictors);
for iNeuron = 1:nNeurons %Loop over all neurons:
    idx_lambda              = output.modelFits_sh(iNeuron).lambda==output.modelFits_sh(iNeuron).(params.lambdastring);
    weightmat(iNeuron,:)    = output.modelFits_sh(iNeuron).glmnet_fit.beta(:,idx_lambda);
end

%normalize the weights:
% maxweights = repmat(nanmax(abs(weightmat),[],1),nNeurons,1); maxweights(maxweights==0) = 1e-6;
% weightmat  = weightmat ./ maxweights;

maxweights = repmat(nanmax(abs(weightmat),[],2),1,params.nPredictors); %maxweights(maxweights==0) = 1e-6;
weightmat  = weightmat ./ maxweights;

g = randperm(nNeurons,nNeurons);
f = randperm(nNeurons,nNeurons);

idx_tri         = ismember(output.x_label,{'trialNum'});
idx_vis         = ismember(output.x_label,{'smallVisualChangepostOri12' 'smallVisualChangepostOri34' 'largeVisualChangepostOri12' 'largeVisualChangepostOri34'});
idx_mot         = ismember(output.x_label,{'videoMotion' 'pupilArea'    'pupilX'    'pupilY'});
idx_aud         = ismember(output.x_label,{'smallAudioChangepostFreq12' 'smallAudioChangepostFreq34' 'largeAudioChangepostFreq12' 'largeAudioChangepostFreq34'});

for iNeuron = 1:nNeurons %Loop over all neurons:
    idx_lambda      = output.modelFits_sh(iNeuron).lambda==output.modelFits_sh(iNeuron).(params.lambdastring);
    output.modelFits(iNeuron).glmnet_fit.beta(:,idx_lambda)             = weightmat(iNeuron,:);
%     output.modelFits(iNeuron).glmnet_fit.beta(idx_tri,idx_lambda)       = 0;
    
%     output.modelFits_sh(iNeuron).glmnet_fit.beta(idx_tri,idx_lambda)    = 0;
    output.modelFits_sh(iNeuron).glmnet_fit.beta(idx_vis,idx_lambda)    = weightmat(iNeuron,idx_vis);
    output.modelFits_sh(iNeuron).glmnet_fit.beta(idx_aud,idx_lambda)    = weightmat(g(iNeuron),idx_aud);
    output.modelFits_sh(iNeuron).glmnet_fit.beta(idx_mot,idx_lambda)    = weightmat(f(iNeuron),idx_mot);
    
%     output.modelFits_sh(iNeuron).glmnet_fit.beta(:,idx_lambda) = weightmat(iNeuron,:);
end

% for iNeuron = 1:nNeurons %Loop over all neurons:
%     idx_lambda = output.modelFits_sh(iNeuron).lambda==output.modelFits_sh(iNeuron).(params.lambdastring);
% %     output.modelFits_sh(iNeuron).glmnet_fit.beta(idx_vis,idx_lambda) = h(iNeuron,idx_vis)*1;
%     output.modelFits_sh(iNeuron).glmnet_fit.beta(idx_aud,idx_lambda) = h(g(iNeuron),idx_aud)*1;
%     output.modelFits_sh(iNeuron).glmnet_fit.beta(idx_mot,idx_lambda) = h(f(iNeuron),idx_mot)*1;
% end

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
    
    if ~isempty(output.modelFits(iNeuron,1).lambda)
        nTotalSpikeBins     = sum(~isnan(output.y(iNeuron,:))); %compute how many total spike bins there are from output struct
        nTrials             = nTotalSpikeBins / params.nTimebins; %how many trials iin this session
        idx_ses             = strcmp(sessionData.session_ID,spikeData.session_ID(iNeuron)); %get which session this is out of all ses predictors
        temptrialData       = MOL_getTempPerSes(sessionData.session_ID(idx_ses),trialData); %filter trialdata for this session
        X_full              = squeeze(output.x(idx_ses,1:nTotalSpikeBins,1:params.nPredictors)); %get predictor matrix for this session
%         Y_full              = squeeze(output.y(iNeuron,1:nTotalSpikeBins))'; %get spike rates
        Yh_full             = cvglmnetPredict(output.modelFits(iNeuron,1),X_full,params.lambdastring,'response'); %cv predicted firing rates for full model
        Yh_full_r           = reshape(Yh_full,params.nTimebins,nTrials)';        
        Yh_sh               = cvglmnetPredict(output.modelFits_sh(iNeuron,1),X_full,params.lambdastring,'response'); %cv predicted firing rates for full model
        Yh_sh_r             = reshape(Yh_sh,params.nTimebins,nTrials)';
    end
    
%     figure;
%     subplot(1,3,1)
%     imagesc(hist_mat)
%     subplot(1,3,2)
%     imagesc(Yh_full_r)
%     subplot(1,3,3)
%     imagesc(Yh_sh_r)

    %Orientation coding: 
    feature     = [ones(sum(splits{1}),1); ones(sum(splits{2}),1)*2]; %Make correct feature vector (ones for condition A and two's for condition B)
    timeidx     = params.xtime>params.twin_resp_start & params.xtime<=params.twin_resp_stop;
    
    responses   = [nanmean(hist_mat(splits{1},timeidx),2); nanmean(hist_mat(splits{2},timeidx),2)]; %get all firing rates in this bin for all trials
    
    if sum(splits{1})>params.min_ntrial && sum(splits{2})>params.min_ntrial %for both conditions the number of trials is met
        [~,~,~,AUC_ORI_V(iNeuron,1)]          = perfcurve2(feature,responses,1); %compute AUC
        AUC_ORI_V(iNeuron,1)                  = AUC_ORI_V(iNeuron,1)*2-1;  %-1 to 1
% 
%         %Make shuffle distribution and compute permutation test threshold:
%         AUC_shuffle = NaN(1,params.nshuffle); %init vector for storing AUC shuffled values
% %         for iShuf = 1:params.nshuffle
%         parfor iShuf = 1:params.nshuffle
%             [~,~,~,AUC_shuffle(iShuf)]    = perfcurve2(feature(randperm(length(feature))),responses,1); %#ok<PFBNS>
%         end
%         
%         AUC_shuffle             = AUC_shuffle*2-1;
%         pVal_ORI_V(iNeuron,1)     = sum(AUC_shuffle>AUC_ORI_V(iNeuron,1))/params.nshuffle;
    end
    
    %Orientation coding during conflict trials: 
    feature     = [ones(sum(splits{3}),1); ones(sum(splits{4}),1)*2]; %Make correct feature vector (ones for condition A and two's for condition B)
    timeidx     = params.xtime>params.twin_resp_start & params.xtime<=params.twin_resp_stop;
    responses   = [nanmean(hist_mat(splits{3},timeidx),2); nanmean(hist_mat(splits{4},timeidx),2)]; %get all firing rates in this bin for all trials
    
    if sum(splits{3})>params.min_ntrial && sum(splits{4})>params.min_ntrial %for both conditions the number of trials is met
        [~,~,~,AUC_ORI_C(iNeuron,1)]          = perfcurve2(feature,responses,1); %compute AUC
        AUC_ORI_C(iNeuron,1)                  = AUC_ORI_C(iNeuron,1)*2-1;  %-1 to 1

        %         %Make shuffle distribution and compute permutation test threshold:
        %         AUC_shuffle = NaN(1,params.nshuffle); %init vector for storing AUC shuffled values
        % %         for iShuf = 1:params.nshuffle
        %         parfor iShuf = 1:params.nshuffle
        %             [~,~,~,AUC_shuffle(iShuf)]    = perfcurve2(feature(randperm(length(feature))),responses,1); %#ok<PFBNS>
        %         end
        %
        %         AUC_shuffle             = AUC_shuffle*2-1;
        %         pVal_ORI_C(iNeuron,1)    = sum(AUC_shuffle>AUC_ORI_C(iNeuron,1))/params.nshuffle;
        
        responses   = [nanmean(Yh_full_r(splits{3},timeidx),2); nanmean(Yh_full_r(splits{4},timeidx),2)]; %get all firing rates in this bin for all trials

        [~,~,~,AUC_ORI_C_MF(iNeuron,1)]          = perfcurve2(feature,responses,1); %compute AUC
        AUC_ORI_C_MF(iNeuron,1)                  = AUC_ORI_C_MF(iNeuron,1)*2-1;  %-1 to 1
        
        responses   = [nanmean(Yh_sh_r(splits{3},timeidx),2); nanmean(Yh_sh_r(splits{4},timeidx),2)]; %get all firing rates in this bin for all trials

        [~,~,~,AUC_ORI_C_MS(iNeuron,1)]          = perfcurve2(feature,responses,1); %compute AUC
        AUC_ORI_C_MS(iNeuron,1)                  = AUC_ORI_C_MS(iNeuron,1)*2-1;  %-1 to 1
        
    end
        
    
    
end

% pVal_ORI_V(pVal_ORI_V==0) = 1/params.nshuffle;
% pVal_ORI_C(pVal_ORI_C==0) = 1/params.nshuffle;


end

