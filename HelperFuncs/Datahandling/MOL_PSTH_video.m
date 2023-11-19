function [hist_mat,hist_mat_tot,hist_mat_z,params] = MOL_PSTH_video(sessionData,trialData,videoData,params)

%Make some dimensional variables:
edges                   = params.t_pre:1e6/25:params.t_post;
params.xtime_video      = edges+1e6/25/2;
params.nTimebins_video  = numel(edges);
nTrials                 = length(trialData.stimChange);

%init variable:
hist_mat                = NaN(params.nTimebins_video,params.nSVDs,nTrials);

fprintf('Computing svd response for trial        \n');
for iTrial = 1:nTrials %for every trial, get the mot svd over time around stim and interpolate to time axis:
    fprintf(repmat('\b', 1, numel([num2str(iTrial-1) num2str(nTrials)])+1));
    fprintf('%d/%d',iTrial,nTrials);
    ses_idx                             = strcmp(sessionData.session_ID,trialData.session_ID(iTrial));
    hist_mat(:,:,iTrial)                = interp1(videoData.ts{ses_idx}-trialData.(params.AlignOn)(iTrial),videoData.(params.videofield){ses_idx}(:,1:params.nSVDs),params.xtime_video,'linear');
end
fprintf('\n')

% Make a z-scored total motion matrix to compare across sessions:
nSessions                       = length(sessionData.session_ID);

hist_mat_tot                    = squeeze(nansum(abs(hist_mat(:,1:params.nSVDs,:)),2));
hist_mat_z                      = hist_mat_tot;

for iSes = 1:nSessions %for every session, subtract baseline and divide by variability during baseline for that session:
    idx_ses                     = strcmp(trialData.session_ID,sessionData.session_ID(iSes));
    temp                        = hist_mat_z(params.xtime_video<0,idx_ses);
    hist_mat_z(:,idx_ses)       = (hist_mat_z(:,idx_ses) - repmat(nanmean(temp,1),params.nTimebins_video,1)) / nanstd(temp(:));
end

end

%Old version: 

% %init variable:
% hist_mat                = NaN(params.nTimebins_video,params.nSVDs,nTrials);
% 
% fprintf('Computing svd response for trial        \n');
% for iTrial = 1:nTrials %for every trial, get the mot svd over time around stim and interpolate to time axis:
%     fprintf(repmat('\b', 1, numel([num2str(iTrial-1) num2str(nTrials)])+1));
%     fprintf('%d/%d',iTrial,nTrials);
%     ses_idx                             = strcmp(sessionData.session_ID,trialData.session_ID(iTrial));
%     hist_mat(:,:,iTrial)                = interp1(videoData.ts{ses_idx}-trialData.(params.AlignOn)(iTrial),videoData.(params.videofield){ses_idx}(:,1:params.nSVDs),params.xtime_video,'linear');
% end
% fprintf('\n')
% 
% % Make a z-scored total motion matrix to compare across sessions:
% nSessions                       = length(sessionData.session_ID);
% 
% hist_mat_tot                    = squeeze(nansum(abs(hist_mat(:,1:params.nSVDs,:)),2));
% hist_mat_z                      = hist_mat_tot;
% 
% for iSes = 1:nSessions %for every session, subtract baseline and divide by variability during baseline for that session:
%     idx_ses                     = strcmp(trialData.session_ID,sessionData.session_ID(iSes));
%     temp                        = hist_mat_z(params.xtime_video<0,idx_ses);
%     hist_mat_z(:,idx_ses)       = (hist_mat_z(:,idx_ses) - repmat(nanmean(temp,1),params.nTimebins_video,1)) / nanstd(temp(:));
% end