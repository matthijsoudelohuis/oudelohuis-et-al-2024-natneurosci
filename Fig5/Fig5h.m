%% Oude Lohuis et al. 2023 Nat Neurosci
% Triple dissociation of auditory, visual, and motor processing in primary visual cortex
% MOL (C) 2023

% Script analyzes the GLM during muscimol and control sessions

startover

%% Parameters for loading: 
%Subselect all MST animals from all fitted GLM sessions: (excluding NE and UST with no muscimol)
params.animals              = {'2003' '2009' '2010' '2011' '2012' '2013' '2019' '2020' '2021' '2022' '2023' '2026' '2027' '2030' '2044' '2045'}; 

%% Load dataset:
% savedate                    = '04-01-23musc';
% folderpath                  = fullfile('E:','Data','Analysis','neuroGLM',savedate);
folderpath                  = fullfile('E:','Matlab','oudelohuis-et-al-2023-natneurosci','Fig5','GLMfits_Musc');
fileList                    = dir(fullfile(folderpath,'*.mat'));
fileList                    = {fileList(:).name};
fileList                    = fileList(~contains(fileList,'X'));
temp                        = cellfun(@(x) x(5:8),fileList,'UniformOutput',false); %get only MST animals
fileList                    = fileList(contains(temp,params.animals));
nFiles                      = length(fileList);

%Load the first session, then concatenate the rest
loadstruct                  = load(fullfile(folderpath,fileList{1}));

output.x_sesid              = cell(nFiles,1);
output.y                    = NaN(2000,size(loadstruct.output.y,2));
output.modelFits            = loadstruct.output.modelFits;
output.x_label              = loadstruct.output.x_label;
nTotalneurons               = 0;

sessionData                 = struct();
trialData                   = struct();
spikeData                   = struct();

for iF = 1:nFiles
    fprintf('Loading and appending data from session #%d/%d\n',iF,nFiles)
    loadstruct          = load(fullfile(folderpath,fileList{iF}));

    sessionData         = AppendStruct(sessionData,loadstruct.sessionData);
    trialData           = AppendStruct(trialData,loadstruct.trialData);
    spikeData           = AppendStruct(spikeData,loadstruct.spikeData);
    
    nNeurons = length(loadstruct.spikeData.session_ID);
    idx = nTotalneurons+1:nTotalneurons+nNeurons;
    
    output.x_sesid(iF,1)= loadstruct.output.x_sesid;
    output.y(idx,:)     = loadstruct.output.y;
    if iF>1
        output.modelFits    = [output.modelFits; loadstruct.output.modelFits];
    end

    nTotalneurons       = nTotalneurons+nNeurons;
   
end

sessionData.Experiment      = strrep(sessionData.Experiment,num2str(2),'');

%% Parameters:
params                      = loadstruct.params;
params                      = MOL_getColors_CHDET(params);
params.savedir              = 'E:\Documents\PhD\Figures\Project CHDET\Results - auV1\7MuscimolA1\GLMMus\';

idx_ses_ctrl                = ~strcmp(sessionData.MuscimolArea,'A1');
idx_ses_musc                = strcmp(sessionData.MuscimolArea,'A1');

idx_ctrl                    = ismember(spikeData.session_ID,sessionData.session_ID(idx_ses_ctrl)) & strcmp(spikeData.area,'V1');
idx_musc                    = ismember(spikeData.session_ID,sessionData.session_ID(idx_ses_musc)) & strcmp(spikeData.area,'V1');
idx_man                     = [idx_ctrl idx_musc];

params.labels_mans          = {'Ctrl' 'Mus'};

%% Report dataset:
fprintf('Dataset: %d sessions, %d trials, %d neurons, %d videos\n',length(sessionData.session_ID),length(trialData.session_ID),length(spikeData.session_ID),length(sessionData.session_ID));
nSessions           = length(sessionData.session_ID);
nNeurons            = length(spikeData.session_ID);
nTrials             = length(trialData.session_ID);

output.y            = output.y(1:nNeurons,:);

%% Multi-level statistics: 
X_coh           = cell(nNeurons,1);
X_coh(ismember(spikeData.session_ID,sessionData.session_ID(ismember(sessionData.Experiment,params.Experiments{1}))))         = params.ExperimentLabels(1);
X_coh(ismember(spikeData.session_ID,sessionData.session_ID(ismember(sessionData.Experiment,params.Experiments{2}))))         = params.ExperimentLabels(2);
X_coh(ismember(spikeData.session_ID,sessionData.session_ID(ismember(sessionData.Experiment,params.Experiments{3}))))         = params.ExperimentLabels(3);

G_mou           = cell(nNeurons,1);
uMice           = unique(sessionData.mousename);
for iMouse = 1:length(uMice)
    G_mou(ismember(spikeData.session_ID,sessionData.session_ID(strcmp(sessionData.mousename,uMice{iMouse})))) = uMice(iMouse);
end

%% Show model prediction for V1 neurons recorded with muscimol in AC

%NOT included in the MS:
close all;
cell_IDs            = {};
cell_IDs{end+1}     = '20312020012421259';
cell_IDs{end+1}     = '20312020012421309';

params.exportfig    = 0;

lastsesid = [];

% for iNeuron = 1:5%params.nNeurons %loop over neurons
for iNeuron = find(ismember(spikeData.cell_ID,cell_IDs))'
    nTotalSpikeBins     = sum(~isnan(output.y(iNeuron,:)));
    nTrials             = nTotalSpikeBins / params.nTimebins;
    
    temptrialData       = MOL_getTempPerSes(spikeData.session_ID(iNeuron),trialData);
    
    if ~(strcmp(spikeData.session_ID(iNeuron),lastsesid))
        load(fullfile(folderpath,sprintf('X_Full_GLM_%s.mat',spikeData.session_ID{iNeuron})))
        lastsesid           = spikeData.session_ID(iNeuron);
    end
    
    Y_full              = squeeze(output.y(iNeuron,1:nTotalSpikeBins))';
    
    Yh_full             = cvglmnetPredict(output.modelFits(iNeuron,1),X_full,params.lambdastring,'response');

    Yvars               = {};
    Yvars{1}            = Y_full;
    Yvars{2}            = Yh_full;
    params.Yvarlabels   = {'True' 'Full model' 'Trial' 'Visual' 'Audio' 'Motor'};
    
    for iC = 1:params.nShuffleCats
        idx             = ismember(output.x_label,params.shuffleVars{iC});
        X_temp          = X_full; 
        X_temp(:,~idx) = 0;  
        Yh              = cvglmnetPredict(output.modelFits(iNeuron,1),X_temp,params.lambdastring,'response');
        Yvars{2+iC}     = Yh;
    end
    MOL_plotYYhat(temptrialData,params,Yvars,spikeData.cell_ID{iNeuron})
end

%% Compute variance explained over all individual trials:
idx_time            = params.xtime>0 & params.xtime<=0.2e6;

cvR2_full           = NaN(nNeurons,1);
cvR2_cat            = NaN(nNeurons,params.nShuffleCats);

fprintf('Computing cvR2 for neuron        \n');

lastsesid = [];

for iNeuron = 1:nNeurons %loop over neurons
    fprintf(repmat('\b', 1, numel([num2str(iNeuron-1) num2str(nNeurons)])+2));
    fprintf('%d/%d\n',iNeuron,nNeurons);
    if ~isempty(output.modelFits(iNeuron,1).lambda)
        nTotalSpikeBins     = sum(~isnan(output.y(iNeuron,:))); %compute how many total spike bins there are from output struct
        nTrials             = nTotalSpikeBins / params.nTimebins; %how many trials iin this session
        temptrialData       = MOL_getTempPerSes(spikeData.session_ID(iNeuron),trialData);
        
        if ~(strcmp(spikeData.session_ID(iNeuron),lastsesid))
            load(fullfile(folderpath,sprintf('X_Full_GLM_%s.mat',spikeData.session_ID{iNeuron})))
            lastsesid           = spikeData.session_ID(iNeuron);
        end
        
        Y_full              = squeeze(output.y(iNeuron,1:nTotalSpikeBins))'; %get spike rates
        Yh_full             = cvglmnetPredict(output.modelFits(iNeuron,1),X_full,params.lambdastring,'response'); %cv predicted firing rates for full model
        idx_time_all        = repmat(idx_time,1,nTrials); %make an index of all time bins that are included for computation of cvR2
        cvR2_full(iNeuron,1) = 1 - var(Y_full(idx_time_all) - Yh_full(idx_time_all)) / var(Y_full(idx_time_all)); %compute R2

        for iC = 1:params.nShuffleCats
            %Predict firing rate when using only one predictor set, compute R2 and store value:
            idx                     = ismember(output.x_label,params.shuffleVars{iC});
            X_temp                  = X_full;
            X_temp(:,~idx)          = 0;
            Yh_cat                  = cvglmnetPredict(output.modelFits(iNeuron,1),X_temp,params.lambdastring,'response');
            cvR2_cat(iNeuron,iC)    = 1 - var(Y_full(idx_time_all) - Yh_cat(idx_time_all)) / var(Y_full(idx_time_all));
        end
    end
end

%% Correct for negative cvR2
cvR2_full(cvR2_full<0)                  = 0;
cvR2_cat(cvR2_cat<0)                    = 0;

%% Show cvR2 for each predictor subset for control and muscimol sessions:
params.labels_cats      = {'Trial' 'Vis' 'Aud' 'Motor'};

figure; set(gcf,'units','normalized','Position',[0.2 0.5 0.2 0.2],'color','w'); hold all;

params.colors_splits = {[0.4 0.4 0.4] [0 0 0.8] [0.8 0 0] [0.1 0.7 0.3]};
handles = [];

fprintf('Variance explained between ctrl and mus sessions: \n')

for iC = 1:params.nShuffleCats
    for iMan = 1:2
        idx     = idx_man(:,iMan);
        tmp     = cvR2_cat(idx,iC);
%         tmp(tmp<0) = 0;

        xloc = iC/2 + iMan/5;
        if iMan == 1
            handles(iC) = bar(xloc,nanmean(tmp),0.18,'k');
            set(handles(iC),'FaceColor',params.colors_splits{iC})
        elseif iMan==2
            handles2(iC) = bar(xloc,nanmean(tmp),0.18,'k');
            set(handles2(iC),'FaceColor',mean([params.colors_splits{iC}; 0.8 0.8 0.8;  0.8 0.8 0.8],1))
        end

        errorbar(xloc,nanmean(tmp),nanstd(tmp)/sqrt(sum(idx)),'k','LineWidth',1,'CapSize',6);
    end
    
    Y                   = cvR2_cat(:,iC);
    xlocs               = iC/2 + [1 2]/5;

    %statistics:
    X_man               = zeros(size(Y));
    X_man(idx_ctrl)     = 1;
    X_man(idx_musc)     = 2;
    
    idx                 = idx_ctrl | idx_musc;
    %statistics:
    fprintf('%s: \n',params.labels_cats{iC})
    tbl             = table(Y(idx),X_man(idx),G_mou(idx),'VariableNames',{'Y','X_man','Mouse'}); %Create table for mixed model
    lme             = fitlme(tbl,'Y~X_man'); %construct linear mixed effects model with fixed effect of temporal window
    stats           = dataset2table(anova(lme)); %Perform ANOVA on model and output as matrix
    
    fprintf('F(%d,%2.0f)=%1.2f, p=%1.3f; \n',stats{2,3},stats{2,4},stats{2,2},stats{2,5})
    sigstar(xlocs,stats{2,5})

end
set(gca,'XTick',1.5:1:5.5,'XTickLabels',params.labels_mans,'YTick',0:0.02:0.2)
xlim([0.5 2.5])
ylim([0 0.07])
legend(handles(1:4),{'Trial' 'Visual' 'Audio' 'Motor'},'Location','best'); legend boxoff;
export_fig(fullfile(params.savedir,sprintf('Bar_cvR2_V1_CtrlMus')),'-eps','-nocrop')

%% Show cvR2 for each predictor subset for control and muscimol sessions:
params.labels_cats      = {'Trial' 'Vis' 'Aud' 'Motor'};

figure; set(gcf,'units','normalized','Position',[0.2 0.5 0.17 0.24],'color','w'); hold all;

params.colors_splits = {[0.4 0.4 0.4] [0 0 0.8] [0.8 0 0] [0.1 0.7 0.3]};
handles = [];

fprintf('Variance explained between ctrl and mus sessions: \n')

datamat = NaN(nNeurons,params.nShuffleCats*2);
% datamat = NaN(params.nShuffleCats*2,nNeurons);
for iC = 1:params.nShuffleCats
    for iMan = 1:2
        idx     = idx_man(:,iMan);
        tmp     = cvR2_cat(idx,iC);
        tmp = tmp(tmp>prctile(tmp,2.5) & tmp<prctile(tmp,95));
        datamat(1:length(tmp),(iC-1)*2+iMan) = tmp;
    end
end

h = violinplot(datamat,[params.labels_cats params.labels_cats],'EdgeColor',[1 1 1],'BandWidth',0.01,...
    'Width',0.3,'ShowData',false,'ViolinAlpha',1);

for iC = 1:params.nShuffleCats
    for iMan = 1:2
        if iMan == 1
            h((iC-1)*2+iMan).ViolinColor = params.colors_splits{iC};
        else
            h((iC-1)*2+iMan).ViolinColor  = mean([params.colors_splits{iC}; 0.8 0.8 0.8;  0.8 0.8 0.8],1);
        end
    end
end

for iC = 1:params.nShuffleCats

    Y                   = cvR2_cat(:,iC);
    
    xlocs               = [1 2] + (iC-1)*2;

    %statistics:
    X_man               = zeros(size(Y));
    X_man(idx_ctrl)     = 1;
    X_man(idx_musc)     = 2;
    
    idx                 = idx_ctrl | idx_musc;
    %statistics:
    fprintf('%s: \n',params.labels_cats{iC})
    tbl             = table(Y(idx),X_man(idx),G_mou(idx),'VariableNames',{'Y','X_man','Mouse'}); %Create table for mixed model
    lme             = fitlme(tbl,'Y~X_man'); %construct linear mixed effects model with fixed effect of temporal window
    stats           = dataset2table(anova(lme)); %Perform ANOVA on model and output as matrix
    
    fprintf('F(%d,%2.0f)=%1.2f, p=%1.3f; \n',stats{2,3},stats{2,4},stats{2,2},stats{2,5})
    sigstar(xlocs,stats{2,5})

end

set(gca,'XTick',1:8,'XTickLabels',params.labels_mans,'YTick',0:0.05:0.2)

xlim([0.5 8.5])
ylim([-0.005 0.2])

export_fig(fullfile(params.savedir,sprintf('Violin_cvR2_V1_CtrlMus')),'-eps','-nocrop')

tbl             = table(cvR2_cat(idx,1),cvR2_cat(idx,2),cvR2_cat(idx,3),cvR2_cat(idx,4),X_man(idx),G_mou(idx),'VariableNames',...
    {'Trial','Vis','Aud','Mot','Manipulation','Mouse'}); %Create table for mixed model
writetable(tbl,'SourceData_Fig5h_EV_Muscimol.xlsx')









