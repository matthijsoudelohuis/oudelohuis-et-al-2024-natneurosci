
%%
startover

%% Load dataset:
% savedate                    = '14-01-22';
savedate                    = '03-02-22';
folderpath                  = fullfile('E:','Data','Analysis','neuroGLM',savedate);
fileList                    = dir(fullfile(folderpath,'*.mat'));
fileList                    = {fileList(:).name};
fileList                    = fileList(cellfun(@isempty,strfind(fileList,'X')));
nFiles                      = length(fileList);

%Load the first session, then concatenate the rest
loadstruct          = load(fullfile(folderpath,fileList{1}));

% output.x            = NaN(nFiles,size(loadstruct.output.x,2),size(loadstruct.output.x,3));
output.x_sesid      = cell(nFiles,1);
output.y            = NaN(2000,size(loadstruct.output.y,2));
output.modelFits    = loadstruct.output.modelFits;
output.shuffleFits  = loadstruct.output.shuffleFits;
output.x_label      = loadstruct.output.x_label;
nTotalneurons       = 0;

sessionData_GLM     = struct();
trialData_GLM       = struct();
spikeData_GLM       = struct();
% videoData           = struct();

for iF = 1:nFiles
    fprintf('Loading and appending data from session #%d/%d\n',iF,nFiles)
    loadstruct          = load(fullfile(folderpath,fileList{iF}));

    sessionData_GLM         = AppendStruct(sessionData_GLM,loadstruct.sessionData);
    trialData_GLM           = AppendStruct(trialData_GLM,loadstruct.trialData);
    spikeData_GLM           = AppendStruct(spikeData_GLM,loadstruct.spikeData);
%     videoData           = AppendStruct(videoData,loadstruct.videoData);
    
    nNeurons = length(loadstruct.spikeData.session_ID);
    idx = nTotalneurons+1:nTotalneurons+nNeurons;
    
%     output.x(iF,:,:)    = loadstruct.output.x;
    output.x_sesid(iF,1)= loadstruct.output.x_sesid;
    output.y(idx,:)     = loadstruct.output.y;
    if iF>1
        output.modelFits    = [output.modelFits; loadstruct.output.modelFits];
        output.shuffleFits  = [output.shuffleFits; loadstruct.output.shuffleFits];
    end
    
    nTotalneurons       = nTotalneurons+nNeurons;
end

%% Parameters:
params                      = loadstruct.params;
params                      = MOL_getColors_CHDET(params);
params.savedir              = 'E:\Documents\PhD\Figures\Project CHDET\Results - auV1\18GLM\Conflicts\';
clear loadstruct;

%%

%% Get data:
% [Data] = MOL_GetData('E:','CHDET',params.Experiments,{'2012' '2031'},{},{'sessionData' 'trialData_newtrials' 'spikeData' 'videoData'});
[Data] = MOL_GetData('E:','CHDET',params.Experiments,{},{},{'sessionData' 'trialData_newtrials' 'spikeData' 'videoData'});
sessionData     = Data.sessionData;
trialData       = Data.trialData;
spikeData       = Data.spikeData;
videoData       = Data.videoData;

%% Remove last n trials:
trialData = MOL_RemoveLastnTrials(trialData,5);

%% Filter only sessions with fitted GLM:
[sessionData,trialData,spikeData,videoData]       = MOL_getTempPerSes(sessionData_GLM.session_ID,sessionData,trialData,spikeData,videoData);

%% Filter only neurons from GLM dataset:
spikeFields     = fieldnames(spikeData);
idx             = ismember(spikeData.cell_ID,spikeData_GLM.cell_ID);
fprintf('Filtered %d/%d neurons from GLM\n',sum(idx),length(spikeData.session_ID));
for iF = 1:length(spikeFields)
    spikeData.(spikeFields{iF}) = spikeData.(spikeFields{iF})(idx);
end

%% Filter neurons on area:
spikeFields     = fieldnames(spikeData);
idx             = ismember(spikeData.area,'V1');
fprintf('Filtered %d/%d neurons based on area\n',sum(idx),length(spikeData.session_ID));
for iF = 1:length(spikeFields)
    spikeData.(spikeFields{iF}) = spikeData.(spikeFields{iF})(idx);
end

%% Filter out trials with photostimulation in V1:
sesids              = sessionData.session_ID(sessionData.UseOpto & strcmp(sessionData.PhotostimArea,'V1'));
idx                 = trialData.hasphotostim==1 & ismember(trialData.session_ID,sesids);
trialFields         = fieldnames(trialData);
fprintf('Removed %d/%d trials  with photostim in V1\n',sum(idx),length(trialData.session_ID));
for iF = 1:length(trialFields)
    trialData.(trialFields{iF}) = trialData.(trialFields{iF})(~idx);
end

%% Filter only sessions with conflict trials:
sesids = unique(trialData.session_ID(strcmp(trialData.trialType,'C')));
[sessionData,trialData,spikeData,videoData]       = MOL_getTempPerSes(sesids,sessionData,trialData,spikeData,videoData);

%% Report dataset:
nSessions           = length(sessionData.session_ID);
nNeurons            = length(spikeData.session_ID);
nTrials             = length(trialData.session_ID);
nVideos             = length(videoData.session_ID);

output.y            = output.y(1:nNeurons,:);
fprintf('Dataset: %d sessions, %d trials, %d V1 neurons, %d videos\n',nSessions,nTrials,nNeurons,nVideos);

for iExp = 1:3
    fprintf('%s: %d\n',params.ExperimentLabels{iExp},sum(strcmp(sessionData.Experiment,params.Experiments(iExp))));
end

%% 





%% Show model prediction of conflict trials:
close all;
cell_IDs            = {};

%visually driven example cell:
cell_IDs{end+1}     = '20302020011621436';
cell_IDs{end+1}     = '20442021042411333';
cell_IDs{end+1}     = '20302020011621288';
cell_IDs{end+1}     = '20442021042811300';
% 
cell_IDs{end+1}     = '20302020011621092';
cell_IDs{end+1}     = '20112018081011125';

cell_IDs{end+1}     = '20442021042811376';
cell_IDs{end+1}     = '20302020011621372';

% cell_IDs{end+1}     = '20442021042811372';
% cell_IDs{end+1}     = '20112018081011151';
% cell_IDs{end+1}     = '20302020011621369';
% cell_IDs{end+1}     = '20312020012221145';
% cell_IDs{end+1}     = '20112018081011106'; %perhaps
% cell_IDs{end+1}     = '20112018081011160'; %perhaps
% cell_IDs{end+1}     = '20302020011621289';
% cell_IDs{end+1}     = '20222019062711254';

% cell_IDs{end+1}     = '10092019030831134';
% cell_IDs{end+1}     = '20122018081431146';

% cell_IDs{end+1}     = '20122018081311186';%perhaps
cell_IDs{end+1}     = '20122018081431146';
% cell_IDs{end+1}     = '20122018081431167'; %au only
% cell_IDs{end+1}     = '20122018081431176'; %perpaps but low firing

% cell_IDs = spikeData.cell_ID(3:4);

params.exportfig    = 0;

lastsesid = [];

cell_IDs     = {'20122018081431146'};

% for iNeuron = 1:5%params.nNeurons %loop over neurons
for iNeuron = find(ismember(spikeData.cell_ID,cell_IDs))'
    if ~(strcmp(spikeData.session_ID(iNeuron),lastsesid))
        %Get the relevant data for each session individually:
        [tempsessionData, temptrialData, tempvideoData]      = MOL_getTempPerSes(spikeData.session_ID(iNeuron),sessionData,trialData,videoData);
        
        [X_full,params] = make_X_full_Conflicts(params,tempsessionData, temptrialData, tempvideoData);
        
        lastsesid           = spikeData.session_ID(iNeuron);
    end
    
    % To be predicted spike count bins:
    events_ts           = temptrialData.(params.AlignOn);
    hist_mat            = calc_psth(events_ts,spikeData.ts{iNeuron},params);    %Construct histogram matrix:
    
    Y_full              = reshape(hist_mat(:,:)',params.nTotalTimebins,1); %linearized vector of all spike count bins
    
    idx_neuron_GLM      = strcmp(spikeData_GLM.cell_ID,spikeData.cell_ID{iNeuron});
    
    Yh_full             = cvglmnetPredict(output.modelFits(idx_neuron_GLM,1),X_full,params.lambdastring,'response');
    
    Yvars               = {};
    
    Yvars{1}            = Y_full;
    Yvars{2}            = Yh_full;
    params.Yvarlabels   = {'True' 'Full model' 'Trial' 'Visual' 'Audio' 'Motor'};
    
    for iC = 1:params.nShuffleCats
        idx             = ismember(output.x_label,params.shuffleVars{iC});
        X_temp          = X_full; 
        X_temp(:,~idx) = 0;  
        Yh              = cvglmnetPredict(output.modelFits(idx_neuron_GLM,1),X_temp,params.lambdastring,'response');
        Yvars{2+iC}     = Yh;
    end
    MOL_plotYYhat_C(temptrialData,params,Yvars,spikeData.cell_ID{iNeuron})
end


%% Compute variance explained over all individual trials:
idx_time            = params.xtime>0 & params.xtime<=0.5e6;

nSplits             = 4;
cvR2_full           = NaN(nNeurons,nSplits);
cvR2_cat            = NaN(nNeurons,nSplits,params.nShuffleCats);

fprintf('Computing cvR2 for neuron        \n');

lastsesid = [];

for iNeuron = 1:nNeurons %loop over neurons
    fprintf(repmat('\b', 1, numel([num2str(iNeuron-1) num2str(nNeurons)])+2));
    fprintf('%d/%d\n',iNeuron,nNeurons);
    
    idx_neuron_GLM      = strcmp(spikeData_GLM.cell_ID,spikeData.cell_ID{iNeuron});
    
    temp = find(idx_neuron_GLM);
    if numel(temp)>1
        idx     = strcmp(spikeData_GLM.session_ID(idx_neuron_GLM),spikeData.session_ID{iNeuron});
        idx_neuron_GLM(temp(~idx)) = false;
    end
    
    if ~isempty(output.modelFits(idx_neuron_GLM,1).lambda)
        
        if ~(strcmp(spikeData.session_ID(iNeuron),lastsesid))
            %Get the relevant data for each session individually:
            [tempsessionData, temptrialData, tempvideoData]      = MOL_getTempPerSes(spikeData.session_ID(iNeuron),sessionData,trialData,videoData);
            
            [X_full,params] = make_X_full_Conflicts(params,tempsessionData, temptrialData, tempvideoData);
            
            lastsesid           = spikeData.session_ID(iNeuron);
        end
        
        % To be predicted spike count bins:
        events_ts           = temptrialData.(params.AlignOn);
        hist_mat            = calc_psth(events_ts,spikeData.ts{iNeuron},params);    %Construct histogram matrix:
        
        Y_full              = reshape(hist_mat(:,:)',params.nTotalTimebins,1); %linearized vector of all spike count bins
        
        Yh_full             = cvglmnetPredict(output.modelFits(idx_neuron_GLM,1),X_full,params.lambdastring,'response');
        
        splits              = {};
        splits{1}           = strcmp(temptrialData.trialType,'P');
        splits{2}           = strcmp(temptrialData.trialType,'X');
        splits{3}           = strcmp(temptrialData.trialType,'Y');
        splits{4}           = strcmp(temptrialData.trialType,'C');
        
        for iSplit = 1:nSplits
            idx_trials                  = repmat(splits{iSplit},1,params.nTimebins);
            idx_trials                  = reshape(idx_trials',params.nTotalTimebins,1);
            idx_time_all                = repmat(idx_time,1,params.nTrials)'; 
            idx_all                     = idx_trials & idx_time_all; %combine timewindow and trials: make an index of all time bins that are included for computation of cvR2

%             cvR2_full(iNeuron,iSplit)   = var(Yh_full(idx_all))/var(Y_full(idx_all)); %compute R2
            cvR2_full(iNeuron,iSplit)   = 1 - var(Y_full(idx_all) - Yh_full(idx_all)) / var(Y_full(idx_all)); %compute R2

            for iC = 1:params.nShuffleCats
                %Predict firing rate when using only one predictor set, compute R2 and store value:
                idx                     = ismember(output.x_label,params.shuffleVars{iC});
                X_temp                  = X_full;
                X_temp(:,~idx)          = 0;
                Yh_cat                  = cvglmnetPredict(output.modelFits(idx_neuron_GLM,1),X_temp,params.lambdastring,'response');
%                 cvR2_cat(iNeuron,iSplit,iC)    = var(Yh_cat(idx_all))/var(Y_full(idx_all));
                cvR2_cat(iNeuron,iSplit,iC)    = 1 - var(Y_full(idx_all) - Yh_cat(idx_all)) / var(Y_full(idx_all));
            end
        end
    end
end


%% Show cvR2 for each predictor subset and area:
params.nExperiments     = length(params.Experiments);

params.nSplits          = 4;
params.colors_splits    = {[0.4 0.4 0.4] [0 0 0.8] [0.8 0 0] [0.1 0.7 0.3]};

params                  = MOL_getColors_CHDET(params);

params.labels_splits    = {'C' 'V' 'A' 'AV'};

figure; set(gcf,'units','normalized','Position',[0.2 0.43 0.16 0.3],'color','w'); hold all;

meantoplot = squeeze(nanmean(cvR2_cat(:,:,:),1));

handles = bar(meantoplot,0.7,'stacked');

for i=1:4
    handles(i).FaceColor = params.colors_splits{i};
end
% for i=1:3
%     handles(i).FaceColor = params.colors_splits{i+1};
% end

set(gca,'XTick',1:4,'XTickLabels',params.labels_splits,'YTick',0:0.01:0.2)
xlim([0.5 4.5])
% ylim([0 0.055])
legend(handles,{'Trial' 'Visual' 'Audio' 'Motor'},'Location','NorthWest'); legend boxoff;

export_fig(fullfile(params.savedir,sprintf('StackedBar_cvR2_trialtypes')),'-eps','-nocrop')

%% Show cvR2 for each predictor subset and area for each cohort separately:
params.nExperiments     = length(params.Experiments);

figure; set(gcf,'units','normalized','Position',[0.2 0.43 0.45 0.3],'color','w'); hold all;

params.colors_splits = {[0.4 0.4 0.4] [0 0 0.8] [0.8 0 0] [0.1 0.7 0.3]};
handles = [];

for iExp = 1:3
    subplot(1,3,iExp); hold all;
    idx_exp                 = ismember(spikeData.session_ID,sessionData.session_ID(strcmp(sessionData.Experiment,params.Experiments(iExp))));
    for iSplit = 1:params.nSplits
        idx     = idx_exp;
        
        for iC = 1:params.nShuffleCats
            tmp     = cvR2_cat(idx,iSplit,iC);
            handles(iC) = bar(iC + iSplit/(params.nSplits+1),nanmean(tmp),0.18,'k');
            set(handles(iC),'FaceColor',params.colors_splits{iC})
            errorbar(iC + iSplit/(params.nSplits+1),nanmean(tmp),nanstd(tmp)/sqrt(sum(idx)),'k','LineWidth',1,'CapSize',6);
        end
%         for iC = 1:params.nShuffleCats
%             tmp     = cvR2_cat(idx,iSplit,iC);
%             handles(iC) = bar(iSplit + iC/(params.nShuffleCats+1),nanmean(tmp),0.18,'k');
%             set(handles(iC),'FaceColor',params.colors_splits{iC})
%             errorbar(iSplit + iC/(params.nShuffleCats+1),nanmean(tmp),nanstd(tmp)/sqrt(sum(idx)),'k','LineWidth',1,'CapSize',6);
%         end
    end
    set(gca,'XTick',1.5:1:5.5,'XTickLabels',params.labels_splits,'YTick',0:0.01:0.2)
    xlim([1 1+params.nSplits])
    ylim([0 0.06])
    title(params.ExperimentLabels(iExp),'FontSize',10)
    if iExp == 1
        legend(handles(1:4),{'Trial' 'Visual' 'Audio' 'Motor'},'Location','best'); legend boxoff;
    end
end

%% Show explained variance for each category of predictors over time:

params.labels_cats      = {'Full' 'Trial' 'Vis' 'Aud' 'Motor'};
params.colors_cats      = {[0 0 0] [0.5 0.2 0.3] [0 0 0.8] [0.8 0 0] [0.1 0.7 0.3]};

% iExp = 3;
% idx_exp                 = ismember(spikeData.session_ID,sessionData.session_ID(strcmp(sessionData.Experiment,params.Experiments(iExp))));
idx_exp                 = ismember(spikeData.session_ID,sessionData.session_ID(ismember(sessionData.Experiment,params.Experiments)));

idx_area                = strcmp(spikeData.area,'V1');
idx                     = idx_exp & idx_area;

figure; set(gcf,'units','normalized','Position',[0.05 0.59 0.37 0.25],'color','w'); hold all;
subplot(1,2,1); hold all;
meantoplot      = nanmean(cvR2_full_time(idx,:),1);
errortoplot     = nanstd(cvR2_full_time(idx,:),[],1) / sqrt(sum(idx));
plot(params.xtime,meantoplot,'Color',params.colors_cats{1},'LineWidth',0.5);

for iC = 1:4
    meantoplot      = squeeze(nanmean(cvR2_cat_time(idx,iC,:),1));
%     meantoplot      = squeeze(nanmean(cvR2_cat_un_time(idx,iC,:),1));
    errortoplot     = squeeze(nanstd(cvR2_cat_time(idx,iC,:),[],1) / sqrt(sum(idx)));
    plot(params.xtime,meantoplot,'Color',params.colors_cats{iC+1},'LineWidth',0.5)
end

% meantoplot      = squeeze(nanmean(nansum(cvR2_cat_time(idx,:,:),2),1));
% plot(params.xtime,meantoplot,'Color',[1 0 1],'LineWidth',0.5)

xlim([params.t_pre params.t_post])
ylim([0 0.15])
plot([0 0],get(gca,'ylim'),'k:','LineWidth',0.5)
set(gca,'XTick',-2e6:0.2e6:2e6,'XTickLabels',(-2e6:0.2e6:2e6)*1e-3,'YTick',get(gca,'ylim'));
xlabel('Time (ms)')
title('V1')

idx_area = strcmp(spikeData.area,'A1');
idx         = idx_exp & idx_area;

subplot(1,2,2); hold all;
meantoplot      = nanmean(cvR2_full_time(idx,:),1);
errortoplot     = nanstd(cvR2_full_time(idx,:),[],1) / sqrt(sum(idx));
plot(params.xtime,meantoplot,'Color',params.colors_cats{1},'LineWidth',0.5);

for iC = 1:4
    meantoplot      = squeeze(nanmean(cvR2_cat_time(idx,iC,:),1));
%     meantoplot      = squeeze(nanmean(cvR2_cat_un_time(idx,iC,:),1));
    errortoplot     = squeeze(nanstd(cvR2_cat_time(idx,iC,:),[],1) / sqrt(sum(idx)));
    plot(params.xtime,meantoplot,'Color',params.colors_cats{iC+1},'LineWidth',0.5)
end
xlim([params.t_pre params.t_post])
ylim([0 0.15])
plot([0 0],get(gca,'ylim'),'k:','LineWidth',0.5)
set(gca,'XTick',-2e6:0.2e6:2e6,'XTickLabels',(-2e6:0.2e6:2e6)*1e-3,'YTick',get(gca,'ylim'));
xlabel('Time (ms)')
title('AC')
legend(params.labels_cats,'FontSize',12); legend boxoff;

%% Compute predicted firing rate response for trial types based on different subpredictors:

fprintf('Computing firing rate response based on subpredictors for neuron        \n');

params.nSplits      = 8;
ratemat             = NaN(nNeurons,params.nShuffleCats+2,params.nSplits,params.nTimebins);

%Store for each split whether the 0-200 ms response is different from baseline. Third dimension different predictions:
%Types are: original data, full model, model without movement, original data - movement prediction
signmat             = false(nNeurons,params.nSplits,4); 

params.minTrialCond = 10;

params.twin_baseline_start = -1000e3;
params.twin_baseline_stop = 0;
params.twin_resp_start = 0;
params.twin_resp_stop = 200e3;

params.ttestalpha = 0.025;

for iNeuron = 1:nNeurons %loop over neurons
    fprintf(repmat('\b', 1, numel([num2str(iNeuron-1) num2str(nNeurons)])+2));
    fprintf('%d/%d\n',iNeuron,nNeurons);
    if ~isempty(output.modelFits(iNeuron,1).lambda)
        
        nTotalSpikeBins     = sum(~isnan(output.y(iNeuron,:)));
        nTrials             = nTotalSpikeBins / params.nTimebins;
        
        idx_ses             = strcmp(sessionData.session_ID,spikeData.session_ID(iNeuron));
        temptrialData       = MOL_getTempPerSes(sessionData.session_ID(idx_ses),trialData);
        
        splits              = {};
        splits{1}           = ismember(temptrialData.trialType,{'X'}) & temptrialData.visualOriChangeNorm==2 & temptrialData.visualOriPostChangeNorm==1;
        splits{2}           = ismember(temptrialData.trialType,{'X'}) & temptrialData.visualOriChangeNorm==2 & temptrialData.visualOriPostChangeNorm==2;
        splits{3}           = ismember(temptrialData.trialType,{'X'}) & temptrialData.visualOriChangeNorm==3 & temptrialData.visualOriPostChangeNorm==1;
        splits{4}           = ismember(temptrialData.trialType,{'X'}) & temptrialData.visualOriChangeNorm==3 & temptrialData.visualOriPostChangeNorm==2;
        
        splits{5}           = ismember(temptrialData.trialType,{'Y'}) & temptrialData.audioFreqChangeNorm==2 & temptrialData.audioFreqPostChangeNorm==1;
        splits{6}           = ismember(temptrialData.trialType,{'Y'}) & temptrialData.audioFreqChangeNorm==2 & temptrialData.audioFreqPostChangeNorm==2;
        splits{7}           = ismember(temptrialData.trialType,{'Y'}) & temptrialData.audioFreqChangeNorm==3 & temptrialData.audioFreqPostChangeNorm==1;
        splits{8}           = ismember(temptrialData.trialType,{'Y'}) & temptrialData.audioFreqChangeNorm==3 & temptrialData.audioFreqPostChangeNorm==2;
        
        X_full              = squeeze(output.x(idx_ses,1:nTotalSpikeBins,1:params.nPredictors));
        Y_full              = squeeze(output.y(iNeuron,1:nTotalSpikeBins))';
        
        Y_r                 = reshape(Y_full,params.nTimebins,nTrials);
        
        for iSplit = 1:params.nSplits
            ratemat(iNeuron,1,iSplit,:)     = nanmean(Y_r(:,splits{iSplit}),2);
        end
        
        Yh_full             = cvglmnetPredict(output.modelFits(iNeuron,1),X_full,params.lambdastring,'response');
        Yh_r                = reshape(Yh_full,params.nTimebins,nTrials);
        
        for iSplit = 1:params.nSplits
            ratemat(iNeuron,2,iSplit,:)     = nanmean(Yh_r(:,splits{iSplit}),2);
        end
        
        for iC = 1:params.nShuffleCats
            %Predict firing rate when using only one predictor set, compute R2 and store value:
            idx                     = ismember(output.x_label,params.shuffleVars{iC});
            X_temp                  = X_full;
            X_temp(:,~idx)          = 0;
            Yh_cat                  = cvglmnetPredict(output.modelFits(iNeuron,1),X_temp,params.lambdastring,'response');
            Yh_cat_r                = reshape(Yh_cat,params.nTimebins,nTrials);

            for iSplit = 1:params.nSplits
                ratemat(iNeuron,iC+2,iSplit,:)     = nanmean(Yh_cat_r(:,splits{iSplit}),2);
            end
        end
        
        %Compute fracion of significantly responsive neurons on original as well as when regressing out movement variability:
        %First predict firing rate when using full model expect for motor-predictors:
        iC = 4;
        idx                     = ismember(output.x_label,params.shuffleVars{iC});
        X_temp                  = X_full;
        X_temp(:,idx)           = 0;
        Yh_nomot                = cvglmnetPredict(output.modelFits(iNeuron,1),X_temp,params.lambdastring,'response');
        Yh_nomot_r              = reshape(Yh_nomot,params.nTimebins,nTrials);
        
        X_temp                  = X_full;
        X_temp(:,~idx)          = 0;
        Yh_mot                  = cvglmnetPredict(output.modelFits(iNeuron,1),X_temp,params.lambdastring,'response');
        Yh_mot_r                = reshape(Yh_mot,params.nTimebins,nTrials);

        Yh_submot_r             = Y_r - Yh_mot_r;
        
        %Compute baseline activity using original or predicted data without movement:
        varstrings = {'Y_r' 'Yh_r' 'Yh_nomot_r' 'Yh_submot_r'};
        
        for iVar = 1:size(signmat,3)
            for iSplit = 1:params.nSplits %Store the mean response for each of these splits
                if sum(splits{iSplit})>=params.minTrialCond
                    eval(sprintf('tempdata = %s'';',varstrings{iVar}))
                    bsl                                     = nanmean(tempdata(splits{iSplit},params.xtime>params.twin_baseline_start & params.xtime<params.twin_baseline_stop),2);
                    resp                                    = nanmean(tempdata(splits{iSplit},params.xtime>params.twin_resp_start & params.xtime<params.twin_resp_stop),2);
                    [~,signmat(iNeuron,iSplit,iVar)]        = signrank(bsl,resp,'alpha',params.ttestalpha);
                end
            end
        end
        
    end
end

%% Show mean firing rate response to trial types with predictions based on different subsets. 

params.labels_splits    = {'VthrAB' 'VthrCD' 'VmaxAB' 'VmaxCD' 'AthrAB' 'AthrCD' 'AmaxAB' 'AmaxCD';};
params.labels_cats      = {'Raw' 'Full' 'Trial' 'Vis' 'Aud' 'Motor'};
params.colors_cats      = {[0 0 0] [0.5 0.5 0.5] [0.5 0.2 0.3] [0 0 0.8] [0.8 0 0] [0.1 0.7 0.3]};

% iExp = 3;
%     idx_exp                 = ismember(spikeData.session_ID,sessionData.session_ID(strcmp(sessionData.Experiment,params.Experiments(iExp))));
idx_exp                 = ismember(spikeData.session_ID,sessionData.session_ID(ismember(sessionData.Experiment,params.Experiments)));

idx_area                = strcmp(spikeData.area,'V1');
idx                     = idx_exp & idx_area;

figure; set(gcf,'units','normalized','Position',[0.05 0.19 0.9 0.6],'color','w'); hold all;
for iSplit = 1:params.nSplits
    subplot(2,params.nSplits/2,iSplit); hold all;
    for iC = [1 2 4 5 6]
        
        tmp             = squeeze(nanmean(ratemat(idx,iC,iSplit,:),1));
        
        plot(params.xtime,tmp,'Color',params.colors_cats{iC},'LineWidth',0.5)
    end
    xlim([params.t_pre params.t_post])
    ylim([4 9])
    plot([0 0],get(gca,'ylim'),'k:','LineWidth',0.5)
    set(gca,'XTick',[],'YTick',get(gca,'ylim'));
    if iSplit==params.nSplits
        legend(params.labels_cats([1 2 4 5 6]),'FontSize',8);
        legend boxoff
    end
    title(params.labels_splits{iSplit},'FontSize',9)
end

idx_area = strcmp(spikeData.area,'A1');
idx         = idx_exp & idx_area;

figure; set(gcf,'units','normalized','Position',[0.05 0.19 0.9 0.6],'color','w'); hold all;
for iSplit = 1:params.nSplits
    subplot(2,params.nSplits/2,iSplit); hold all;
    for iC = [1 2 4 5 6]
        
        tmp             = squeeze(nanmean(ratemat(idx,iC,iSplit,:),1));
        
        plot(params.xtime,tmp,'Color',params.colors_cats{iC},'LineWidth',0.5)
    end
    xlim([params.t_pre params.t_post])
    ylim([4 11])
    plot([0 0],get(gca,'ylim'),'k:','LineWidth',0.5)
    set(gca,'XTick',[],'YTick',get(gca,'ylim'));
    if iSplit==params.nSplits
        legend(params.labels_cats([1 2 4 5 6]),'FontSize',8);
        legend boxoff
    end
    title(params.labels_splits{iSplit},'FontSize',9)
end

%%



%% Correlation between cvR2 of auditory, visual, and movement variables:

figure; set(gcf,'units','normalized','Position',[0.05 0.59 0.25 0.17],'color','w'); hold all;
subplot(1,2,1); hold all;

idx_exp                 = ismember(spikeData.session_ID,sessionData.session_ID(ismember(sessionData.Experiment,params.Experiments)));

idx_area                = strcmp(spikeData.area,'V1');
idx                     = idx_exp & idx_area;

xdata                   = cvR2_cat(idx,2);
ydata                   = cvR2_cat(idx,4);

scatter(xdata,ydata,30,'bo','filled','MarkerFaceAlpha',0.4) %Color',params.colors_cats{1},'LineWidth',0.5);

ylim([0 0.4]); xlim([0 0.4])

%statistics:
[r,p] = corr(xdata,ydata,'rows','complete');
fprintf('Variance explained was uncorrelated between visual and motor predictors (r = %1.2f, p = %1.2f), \n',r,p)

subplot(1,2,2); hold all;

xdata                   = cvR2_cat(idx,3);
ydata                   = cvR2_cat(idx,4);

[r,p] = corr(xdata,ydata,'rows','complete');
fprintf(' and weakly correlated between auditory and motor predictors (r = %1.2f, p = %1.3e),\n',r,p)

scatter(xdata,ydata,30,'ro','filled','MarkerFaceAlpha',0.4) %Color',params.colors_cats{1},'LineWidth',0.5);
ylim([0 0.4]); xlim([0 0.1])

export_fig(fullfile(params.savedir,sprintf('Scatter_cvR2_AV_3cohorts')),'-eps','-nocrop')

%%

weightmat = NaN(nNeurons,params.nPredictors+1);

for iNeuron = 1:nNeurons
    if ~isempty(output.modelFits(iNeuron,1).lambda)
        weightmat(iNeuron,:) = glmnetCoef(output.modelFits(iNeuron).glmnet_fit,output.modelFits(iNeuron).(params.lambdastring));
    end
end

weightmat = abs(weightmat);

weightmat_cat = NaN(nNeurons,params.nShuffleCats);

for iC = 1:4
    idx = ismember(output.x_label,params.modelVars{1}(ismember(params.modelVars{1},params.shuffleVars{iC})));
    weightmat_cat(:,iC) = sum(weightmat(:,idx),2);
end


%%
idx_exp                 = ismember(spikeData.session_ID,sessionData.session_ID(ismember(sessionData.Experiment,params.Experiments)));

idx_area                = strcmp(spikeData.area,'V1');
idx                     = idx_exp & idx_area;


temp                    = cvR2_cat(idx,:);

sortval                 = (-(temp(:,2) / max(temp(:,2))) + (temp(:,4) / max(temp(:,4)))) ./ (temp(:,3) / max(temp(:,3)));
[~,sortidx]             = sort(sortval);

% temp(sortidx,:) = temp;
temp = temp(sortidx,:);
figure; set(gcf,'units','normalized','Position',[0.05 0.59 0.16 0.2],'color','w'); hold all;

for iC = 2:4 
    subplot(3,1,iC-1); hold all;
    h = bar(1:sum(idx),temp(:,iC),1.5,'k');
    ylim(round(get(gca,'ylim'),2));
    set(gca,'XTick',[],'YTick',get(gca,'ylim'));
    ylabel('cvR2')
%     title(params.labels_cats(iC+2))
    text(sum(idx)/2.5,0.06,params.labels_cats(iC+2),'FontSize',8)
end

export_fig(fullfile(params.savedir,sprintf('Bars_cvR2_AV_3cohorts')),'-eps','-nocrop')

%%



%% Set as significant if responding to one the two stimuli:

signmat_any(:,1,:) = any(signmat(:,[1 2],:),2);
signmat_any(:,2,:) = any(signmat(:,[3 4],:),2);
signmat_any(:,3,:) = any(signmat(:,[5 6],:),2);
signmat_any(:,4,:) = any(signmat(:,[7 8],:),2);


%% Show fracion of significantly responsive neurons when regressing out movement variability:

params.colors_trialtypes    = [params.colors_visual_opto(1:2) params.colors_audio_opto(1:2)];
params.colors_trialtypes    = [params.colors_visual_opto([2 1]) params.colors_audio_opto([2 1])];

params.colors_trialtypes    = [params.colors_ztrials(1:2) {[0.9 0.3 0.1] [0.5 0 0.1]}];

figure; set(gcf,'units','normalized','Position',[0.05 0.4 0.33 0.33],'color','w'); hold all;

% iExp = 3;
%     idx_exp                 = ismember(spikeData.session_ID,sessionData.session_ID(strcmp(sessionData.Experiment,params.Experiments(iExp))));
idx_exp                 = ismember(spikeData.session_ID,sessionData.session_ID(ismember(sessionData.Experiment,params.Experiments)));


idx_area    = strcmp(spikeData.area,'V1');
idx_area    = strcmp(spikeData.area,'V1') & idx_exp;

idx_models  = [1 4];

labels      =  {'Raw' '-Motor'};

colors_models = {[0 0 0] [0.4 0.6 0.8] [0.3 0.4 0.1] [0.8 0.2 0.6]};

subplot(1,2,1); hold all;
for iSplit = 1:4
    
    datatoplot = squeeze(sum(signmat_any(idx_area,iSplit,idx_models),1)) / sum(idx_area);
    h = bar((1:length(idx_models))+(iSplit-1)*2,datatoplot,'FaceColor',params.colors_trialtypes{iSplit});
end
set(gca,'XTick',1:8,'XTickLabels',repmat(labels,1,4),'YTick',0:0.2:1);
ylim([0 0.6])
title('V1')

idx_area    = strcmp(spikeData.area,'A1');
subplot(1,2,2); hold all;
for iSplit = 1:4
    datatoplot = squeeze(sum(signmat_any(idx_area,iSplit,idx_models),1)) / sum(idx_area);
    h = bar((1:length(idx_models))+(iSplit-1)*2,datatoplot,'FaceColor',params.colors_trialtypes{iSplit});
end
set(gca,'XTick',1:8,'XTickLabels',repmat(labels,1,4),'YTick',0:0.2:1);
ylim([0 0.8])
legend(params.labels_splits)
title('AC')

%% 



%% Depthcorrections:
tempcomb = {'20212019070104' 150; '20232019070170' -100; '20282019121322' -100; '20282019121626' -100; '20292019121240' -50; '20292019121343' -50};

for iSes = 1:size(tempcomb,1)
    spikeData.ChannelY(strcmp(spikeData.session_ID,tempcomb{iSes,1})) = spikeData.ChannelY(strcmp(spikeData.session_ID,tempcomb{iSes,1})) + tempcomb{iSes,2};
end

%%
params.depthcorrection          = -50;

spikeData.ChannelY              = spikeData.ChannelY + params.depthcorrection;

%% Make SUA heatmap over laminar depth for each trial type and for each model split:
%Histogram with binning on depth:

params.chdepthmin       = -12.5;
params.chdepthmax       = 1012.5;
params.resolution       = 25;

% binedges                = 0:25:1025;
binedges                = params.chdepthmin:params.resolution:params.chdepthmax;
nBins                   = length(binedges)-1;

params.finalYaxis       = binedges(1:end-1)+params.resolution/2;

nModelSplits            = size(ratemat,2);
SUAdepthmap_all         = zeros(nModelSplits,params.nSplits,nBins,params.nTimebins);


% iExp = 1;
%     idx_exp                 = ismember(spikeData.session_ID,sessionData.session_ID(strcmp(sessionData.Experiment,params.Experiments(iExp))));
idx_area                = strcmp(spikeData.area,'V1') & idx_exp;
idx_area                = strcmp(spikeData.area,'V1');% & idx_exp;

ratemat_z               = NaN(size(ratemat));
for iNeuron = 1:nNeurons
    bsl = ratemat(iNeuron,:,:,params.xtime<0);
    ratemat_z(iNeuron,:,:,:) = (ratemat(iNeuron,:,:,:) - repmat(mean(bsl,4),1,1,1,params.nTimebins));
    ratemat_z(iNeuron,:,:,:) = ratemat_z(iNeuron,:,:,:) ./ repmat(nanstd(bsl(:)),size(ratemat_z(iNeuron,:,:,:)));
end


for iMod = 1:nModelSplits
    for iSplit = 1:params.nSplits
        for iBin = 1:nBins
            idx = idx_area & spikeData.ChannelY>=binedges(iBin) & spikeData.ChannelY<binedges(iBin+1);
            
            SUAdepthmap_all(iMod,iSplit,iBin,:) = nanmean(ratemat_z(idx,iMod,iSplit,:),1);
        end
    end
end

SUAdepthmap_all(isnan(SUAdepthmap_all)) = 0;

%Filter temporally and spatially:
spat_filter     = fspecial('gaussian',[50 50],2);

for iMod = 1:nModelSplits
    for iSplit = 1:params.nSplits
        temp = squeeze(SUAdepthmap_all(iMod,iSplit,:,:));
        SUAdepthmap_all(iMod,iSplit,:,:) = conv2(temp,spat_filter,'same');
    end
end

%%

plotModels      = [1 4 6 1 5 6];
plotTrialTypes  = {[1 2 3 4] [1 2 3 4] [1 2 3 4] [5 6 7 8] [5 6 7 8] [5 6 7 8]};
% plotTrialTypes  = {[3 4] [3 4] [3 4] [7 8] [7 8] [7 8]};
% plotTrialTypes  = {[1 2] [1 2] [1 2] [3 4] [3 4] [3 4]};
% plotTrialTypes  = {[1 2 3 4] [1 2 3 4] [1 2 3 4] [5 6 7 8] [5 6 7 8] [5 6 7 8]};

figure; set(gcf,'units','normalized','Position',[0.05 0.5 0.4 0.3],'color','w');
timeticks = [-0.2e6 0 0.2e6 0.4e6 0.6e6];
binticks = [200 550 900];

for iMod = 1:length(plotModels)
%     subplot(1,6,iMod);        hold all;
    subplot(2,3,iMod);        hold all;
    
    title(params.labels_cats(plotModels(iMod)),'FontSize',10);
    
    SUAdepthmap = squeeze(nanmean(SUAdepthmap_all(plotModels(iMod),plotTrialTypes{iMod},:,:),2));
    
%     imagesc(flipud(SUAdepthmap));
%     imagesc(params.xtime,params.finalYaxis,flipud(SUAdepthmap)); 
    imagesc(params.xtime,params.finalYaxis,SUAdepthmap); 
    plot([0 0],[-2000 2000],':','Color',[1 1 1],'LineWidth',1)
%     caxis([-0.3 1.6])
%             caxis([-0.3 2])
%     colormap(getPyPlot_cMap('rainbow'))

    %Figure make up:
    if iMod==1 || iMod==4
        ylabel('Cortical depth (\mum)','FontSize', 9)
    end
    if iMod==5
        xlabel('Time (ms)','FontSize',9)
    end
    set(gca,'YDir','reverse')
    set(gca,'XTick',timeticks,'XTickLabel',timeticks*1e-3)
    set(gca,'YTick',binticks,'YTickLabel',binticks,'fontsize',7,'tickdir','out')
    xlim([find(params.xtime>-0.2e6,1) find(params.xtime<0.6e6,1,'last')])
    xlim([-0.2e6 0.6e6]);
    ylim([0 1000])
    
%     if iExp==3 && iMod==2
%         cb=colorbar;
%         cb.Position = cb.Position + 1e-10;
%     end
end

export_fig(fullfile(params.savedir,sprintf('LaminarHeatmap_ModelSplits_AV_3cohorts')),'-eps','-nocrop')







