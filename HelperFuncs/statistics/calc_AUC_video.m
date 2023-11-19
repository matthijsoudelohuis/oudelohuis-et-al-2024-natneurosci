function [AUC_ORI_VID,AUC_FREQ_VID,pVal_ORI_VID,pVal_FREQ_VID] = calc_AUC_video(params,trialData,videoData,video_hist_z)
%% This function calculates the area under the receiver operating characteristic (auROC)
% To identify which variables are encoded in single neurons at which timepoints you can identify how well an
% external observer could discriminate variables from the firing rate at single time points. In other words,
% we used the ROC analysis (Green & Swets, 1966) to identify class-differential responses.
% To identify significant coding, we computed the area under the ROC curve for the firing rate distributions
% between two selections of trials (see below) versus a shuffled distribution.
% Matthijs oude Lohuis 2017-2021

video_hist_z = video_hist_z';

%% Init output fields:
nSessions                   = length(videoData.session_ID);

AUC_ORI_VID                 = NaN(nSessions,1);
AUC_FREQ_VID                = NaN(nSessions,1);

pVal_ORI_VID                = NaN(nSessions,1);
pVal_FREQ_VID               = NaN(nSessions,1);

%% Get the right indices:
splits         = {};
splits{1}      = strcmp(trialData.trialType,'X') & trialData.visualOriChangeNorm==3 & ismember(trialData.visualOriPostChangeNorm,[1 2]);
splits{2}      = strcmp(trialData.trialType,'X') & trialData.visualOriChangeNorm==3 & ismember(trialData.visualOriPostChangeNorm,[3 4]);

splits{3}      = strcmp(trialData.trialType,'Y') & trialData.audioFreqChangeNorm==3 & ismember(trialData.audioFreqPostChangeNorm,[1 2]);
splits{4}      = strcmp(trialData.trialType,'Y') & trialData.audioFreqChangeNorm==3 & ismember(trialData.audioFreqPostChangeNorm,[3 4]);

%% Main loop to get psth matrix:
% lastsesid               = []; %memory var to keep track of whether neuron comes from different session and design matrix needs to be reconstructed

fprintf('Computing AUC for session      \n');
timeidx     = params.xtime_video>params.twin_resp_start & params.xtime_video<=params.twin_resp_stop;

for iSes = 1:nSessions %Loop over all neurons:
    fprintf(repmat('\b', 1, numel([num2str(iSes-1) num2str(nSessions)])+2));
    fprintf('%d/%d\n',iSes,nSessions);
    
    idx_ses = strcmp(trialData.session_ID,videoData.session_ID(iSes));
    
%     if ~strcmp(lastsesid,videoData.session_ID(iSes)) %construct new predictor matrix if neuron comes from a new session:
%         %Get the relevant data for each session individually:
%         temptrialData        = MOL_getTempPerSes(videoData.session_ID(iSes),trialData);
%         lastsesid            = videoData.session_ID(iSes); %save this session_ID
%         
%         %Get the right indices:
%         splits                  = {};
%         splits{1}      = strcmp(temptrialData.trialType,'X') & temptrialData.visualOriChangeNorm==3 & ismember(temptrialData.visualOriPostChangeNorm,[1 2]);
%         splits{2}      = strcmp(temptrialData.trialType,'X') & temptrialData.visualOriChangeNorm==3 & ismember(temptrialData.visualOriPostChangeNorm,[3 4]);
%         
%         splits{3}      = strcmp(temptrialData.trialType,'Y') & temptrialData.audioFreqChangeNorm==3 & ismember(temptrialData.audioFreqPostChangeNorm,[1 2]);
%         splits{4}      = strcmp(temptrialData.trialType,'Y') & temptrialData.audioFreqChangeNorm==3 & ismember(temptrialData.audioFreqPostChangeNorm,[3 4]);
%     end
    
    %Orientation coding from video: 
    feature     = [ones(sum(idx_ses & splits{1}),1); ones(sum(idx_ses & splits{2}),1)*2]; %Make correct feature vector (ones for condition A and two's for condition B)
    responses   = [nanmean(video_hist_z(idx_ses & splits{1},timeidx),2); nanmean(video_hist_z(idx_ses & splits{2},timeidx),2)]; %get all firing rates in this bin for all trials
%     feature     = [ones(sum(splits{1}),1); ones(sum(splits{2}),1)*2]; %Make correct feature vector (ones for condition A and two's for condition B)
%     responses   = [nanmean(video_hist_z(timeidx,idx_ses & splits{1}),1); nanmean(video_hist_z(timeidx,idx_ses & splits{2}),1)]; %get all firing rates in this bin for all trials
    
%     if sum(splits{1})>params.minTrialCond && sum(splits{2})>params.minTrialCond %for both conditions the number of trials is met
    if sum(feature==1)>params.minTrialCond && sum(feature==2)>params.minTrialCond %for both conditions the number of trials is met
        
        %             outputmat_nSpikes(iNeuron,iVar)        = sum(responses~=0);   %store the amount of nonzero responses:
        %             [~,~,~,outputmat(iNeuron,iVar)]        = perfcurve(feature,responses,1); %compute AUC
        [~,~,~,AUC_ORI_VID(iSes,1)]          = perfcurve2(feature,responses,1); %compute AUC
%         AUC_ORI_VID(iSes,1)                  = abs(AUC_ORI_VID(iSes,1)-0.5)  + 0.5; %make absolute
        AUC_ORI_VID(iSes,1)                  = AUC_ORI_VID(iSes,1)*2-1; %-1 to 1
        
        %Make shuffle distribution and compute permutation test threshold:
        AUC_shuffle = NaN(1,params.nshuffle); %init vector for storing AUC shuffled values
%         for iShuf = 1:params.nshuffle
        parfor iShuf = 1:params.nshuffle
            [~,~,~,AUC_shuffle(iShuf)]    = perfcurve2(feature(randperm(length(feature))),responses,1); %#ok<PFBNS>
        end
        
        pVal_ORI_VID(iSes,1) = sum(AUC_shuffle>AUC_ORI_VID(iSes,1))/params.nshuffle;
    end
    
    %Frequency coding from video: 
    feature     = [ones(sum(idx_ses & splits{3}),1); ones(sum(idx_ses & splits{4}),1)*2]; %Make correct feature vector (ones for condition A and two's for condition B)
    responses   = [nanmean(video_hist_z(idx_ses & splits{3},timeidx),2); nanmean(video_hist_z(idx_ses & splits{4},timeidx),2)]; %get all firing rates in this bin for all trials
%     feature     = [ones(sum(splits{3}),1); ones(sum(splits{4}),1)*2]; %Make correct feature vector (ones for condition A and two's for condition B)
%     responses   = [nanmean(hist_mat(splits{3},timeidx),2); nanmean(hist_mat(splits{4},timeidx),2)]; %get all firing rates in this bin for all trials
    
%     if sum(splits{3})>params.minTrialCond && sum(splits{4})>params.minTrialCond %for both conditions the number of trials is met
    if sum(feature==1)>params.minTrialCond && sum(feature==2)>params.minTrialCond %for both conditions the number of trials is met
        %             outputmat_nSpikes(iNeuron,iVar)        = sum(responses~=0);   %store the amount of nonzero responses:
        %             [~,~,~,outputmat(iNeuron,iVar)]        = perfcurve(feature,responses,1); %compute AUC
        [~,~,~,AUC_FREQ_VID(iSes,1)]          = perfcurve2(feature,responses,1); %compute AUC
%         AUC_FREQ_VID(iSes,1)                  = abs(AUC_FREQ_VID(iSes,1)-0.5)  + 0.5; %make absolute
        AUC_FREQ_VID(iSes,1)                  = AUC_FREQ_VID(iSes,1)*2-1; %-1 to 1

        %Make shuffle distribution and compute permutation test threshold:
        AUC_shuffle = NaN(1,params.nshuffle); %init vector for storing AUC shuffled values
%         for iShuf = 1:params.nshuffle
        parfor iShuf = 1:params.nshuffle
            [~,~,~,AUC_shuffle(iShuf)]    = perfcurve2(feature(randperm(length(feature))),responses,1); %#ok<PFBNS>
        end
        
        pVal_FREQ_VID(iSes,1) = sum(AUC_shuffle>AUC_FREQ_VID(iSes,1))/params.nshuffle;
    end
    
end

pVal_ORI_VID(pVal_ORI_VID==0) = 1/params.nshuffle;
pVal_FREQ_VID(pVal_FREQ_VID==0) = 1/params.nshuffle;

end

