%% Oude Lohuis et al. 2023 Nat Neurosci
% Triple dissociation of auditory, visual, and motor processing in primary visual cortex
% MOL (C) 2023

%Script analyzes the results from the GLM fits (performed previously on 
% computing cluster)

startover

%% Load dataset:
folderpath                  = fullfile('E:','Matlab','oudelohuis-et-al-2024-natneurosci-data','GLMfits');
fileList                    = dir(fullfile(folderpath,'*.mat'));
fileList                    = {fileList(:).name};
fileList                    = fileList(~contains(fileList,'X'));
nFiles                      = length(fileList);

%Load the first session, then concatenate the rest
loadstruct          = load(fullfile(folderpath,fileList{1}));

% output.x            = NaN(nFiles,size(loadstruct.output.x,2),size(loadstruct.output.x,3));
output.x_sesid      = cell(nFiles,1);
output.y            = NaN(2000,size(loadstruct.output.y,2));
output.modelFits    = loadstruct.output.modelFits;
% output.shuffleFits  = loadstruct.output.shuffleFits;
output.x_label      = loadstruct.output.x_label;
nTotalneurons       = 0;

sessionData         = struct();
trialData           = struct();
spikeData           = struct();
% videoData           = struct();

cvR2                = [];
cvR2_cat            = [];
cvR2_cross          = [];
cvR2_cross_cat      = [];
    
for iF = 1:nFiles
    fprintf('Loading and appending data from session #%d/%d\n',iF,nFiles)
    loadstruct          = load(fullfile(folderpath,fileList{iF}));

    sessionData         = AppendStruct(sessionData,loadstruct.sessionData);
    trialData           = AppendStruct(trialData,loadstruct.trialData);
    spikeData           = AppendStruct(spikeData,loadstruct.spikeData);
%     videoData           = AppendStruct(videoData,loadstruct.videoData);
    
    nNeurons = length(loadstruct.spikeData.session_ID);
    idx = nTotalneurons+1:nTotalneurons+nNeurons;
    
    output.x_sesid(iF,1)= loadstruct.output.x_sesid;
    output.y(idx,:)     = loadstruct.output.y;
    if iF>1
        output.modelFits    = [output.modelFits; loadstruct.output.modelFits];
    end
    
    cvR2            = [cvR2; loadstruct.cvR2]; %#ok<*AGROW>
    cvR2_cat        = [cvR2_cat; loadstruct.cvR2_cat];
    cvR2_cross      = [cvR2_cross; loadstruct.cvR2_cross];
    cvR2_cross_cat  = [cvR2_cross_cat; loadstruct.cvR2_cross_cat];
    
    nTotalneurons       = nTotalneurons+nNeurons;
end

sessionData.Experiment      = strrep(sessionData.Experiment,num2str(2),'');

%% Parameters:
params                      = loadstruct.params;
params                      = MOL_getColors_CHDET(params);
params.savedir              = 'E:\OneDrive\PhD\Figures\Project CHDET\Results - auV1\18GLM\';
params.Yvarlabels           = {'True' 'Full model' 'Trial' 'Visual' 'Audio' 'Motor'};

%% %remove animal 1012 with recordings in LM according to histology
spikeData.area(strcmp(spikeData.session_ID,'10122019041022') & strcmp(spikeData.area,'A1')) = deal({'LM'});
spikeData.area(strcmp(spikeData.session_ID,'10122019041126') & strcmp(spikeData.area,'A1')) = deal({'LM'});

%% Report dataset:
nSessions           = length(sessionData.session_ID);
nNeurons            = length(spikeData.session_ID);
nTrials             = length(trialData.session_ID);

output.y            = output.y(1:nNeurons,:);
fprintf('Dataset: %d sessions, %d trials, %d neurons, %d videos\n',nSessions,nTrials,nNeurons,length(sessionData.session_ID));

for iExp = 1:3
    fprintf('%s: %d\n',params.ExperimentLabels{iExp},sum(strcmp(sessionData.Experiment,params.Experiments(iExp))));
end
fprintf('%d V1 and %d AC neurons \n',sum(strcmp(spikeData.area,'V1')),sum(strcmp(spikeData.area,'A1')));

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

%% FIGURES: 

%% Inspiration for Figure 3c: show schematic of model:
figure; set(gcf,'units','normalized','Position',[0.05 0.59 0.16 0.2],'color','w'); hold all;
counter     = 1;
piece       = 5450:5650;

exNeuron        =  '20442021042411333';
exNeuron_idx    =  find(strcmp(spikeData.cell_ID,'20442021042411333'));
exNeuron_idx    =  exNeuron_idx(1);

load(fullfile(folderpath,sprintf('X_Full_GLM_%s.mat',spikeData.session_ID{exNeuron_idx})))

idx = [1 44:-1:11 398+24:-1:398 202+24:-1:202]; %new model

colors = ['k' repmat('g',1,34) repmat('r',1,25) repmat('b',1,25)];

for i = idx
    temp = X_full(piece,i); temp = temp-min(temp); temp = temp / max(temp);
    temp = temp*0.8+0.1;
    plot(temp+counter,colors(counter),'LineWidth',0.2);
    counter = counter+1;
end

temp = smooth(output.y(exNeuron_idx,piece),3); temp = temp-min(temp); temp = temp / max(temp);
temp = temp*3+0.1;
plot(temp+counter,'k','LineWidth',0.5);
counter = counter+1;

set(gca,'XColor','none','YColor','none')

%% Extended data fig 5: Show model prediction:
close all;
cell_IDs            = {};

%Auditory + motor driven example cell:
cell_IDs{end+1}     = '20122018081431146';
%Auditory no motor example cell:
cell_IDs{end+1}     = '10122019041031329'; 
%Motor no aud example cell:
cell_IDs{end+1}     = '20122018081311186';

params.exportfig    = 0;

lastsesid = [];

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
    
    for iC = 1:params.nShuffleCats
        idx             = ismember(output.x_label,params.shuffleVars{iC});
        X_temp          = X_full; 
        X_temp(:,~idx)  = 0;  
        Yh              = cvglmnetPredict(output.modelFits(iNeuron,1),X_temp,params.lambdastring,'response');
        Yvars{2+iC}     = Yh;
    end
    
    MOL_plotYYhat(temptrialData,params,Yvars,spikeData.cell_ID{iNeuron})
end

%% Correct for negative R2: 
cvR2_nn                      = cvR2;
cvR2_nn(cvR2_nn<0)           = 0;
cvR2_cat_nn                  = cvR2_cat;
cvR2_cat_nn(cvR2_cat_nn<0)   = 0;

%% Report variance explained:
idx     = strcmp(spikeData.area,'V1');
fprintf('\nVariance explained across all V1 neurons (%2.3f, IQR %2.3f-%2.3f) (on trials excluding conflict trials) \n',...
    nanmedian(cvR2_nn(idx,1)),prctile(cvR2_nn(idx,1),[25 75]));


%% Show cvR2 for each predictor subset and area: 
params.nExperiments     = length(params.Experiments);
params.areas            = {'V1' 'A1'};% 'PPC' 'CG1'};
params.nAreas           = length(params.areas); 

figure; set(gcf,'units','normalized','Position',[0.2 0.43 0.22 0.3],'color','w'); hold all;

params.colors_splits = {[0.4 0.4 0.4] [0 0 0.8] [0.8 0 0] [0.1 0.7 0.3]};
handles = [];

for iArea = 1:params.nAreas
    idx     = strcmp(spikeData.area,params.areas{iArea});
    for iC = 1:params.nShuffleCats
        tmp     = cvR2_cat_nn(idx,iC);
        
%         (boxplot)
%         tmp     = tmp(tmp~=0);
        h = boxplot(tmp, 'plotstyle','compact','positions',iArea + iC/(params.nShuffleCats+1),...
        'medianstyle','line','boxstyle','outline','outliersize',0.01,'whisker',1,...
        'colors','k','widths',0.18);
        handles(iC) = h(5);
        %         (barplot)
%         handles(iC) = bar(iArea + iC/(params.nShuffleCats+1),nanmean(tmp),0.18,'k');
%         set(handles(iC),'FaceColor',params.colors_splits{iC})
%         errorbar(iArea + iC/(params.nShuffleCats+1),nanmean(tmp),nanstd(tmp)/sqrt(sum(idx)),'k','LineWidth',1,'CapSize',6);
    end
end
h = findobj(gca,'tag','Outliers');
delete(h)
h = findobj(gca,'Tag','Box');
for j=1:length(h)
    patch(get(h(j),'XData'),get(h(j),'YData'),params.colors_splits{mod(j-1,4)+1},'FaceAlpha',1);
end

%statistics:
Y               = cvR2_cat_nn(:);
X_var           = repmat(params.Yvarlabels(3:6),nNeurons,1);
X_var           = X_var(:);
X_area          = repmat(spikeData.area,1,4);
X_area          = X_area(:);
G_mou2          = repmat(G_mou,1,4);
G_mou2          = G_mou2(:);

idx_V1          = strcmp(X_area,'V1');
idx_AC          = strcmp(X_area,'A1');

tbl             = table(Y(idx_V1),X_var(idx_V1),G_mou2(idx_V1),'VariableNames',{'Y','Var','Mouse'}); %Create table for mixed model
lme             = fitlme(tbl,'Y~Var+(1|Mouse)'); %construct linear mixed effects model with fixed effect of temporal window (NO random intercept for mice as cohorts are tested)
stats           = dataset2table(anova(lme,'DFMethod','Satterthwaite')); %Perform ANOVA on model and output as matrix
fprintf('F(%d,%2.0f)=%1.2f, p=%1.2e; \n',stats{2,3},stats{2,4},stats{2,2},stats{2,5})

contrasts       = {[0 1 0 0] [0 0 1 0] [0 0 0 1] [0 1 -1 0] [0 1 0 -1] [0 0 1 -1]};
positions       = {[1.2 1.4] [1.2 1.6] [1.2 1.8] [1.4 1.6] [1.4 1.8] [1.6 1.8]};
contrastlabels  = {'T vs. V' 'T vs. A' 'T vs. M' 'V vs. A' 'V vs. M' 'A vs. M'};

fprintf('Posthoc comparison:\n')
for iC = 1:length(contrasts)
    [p,F,DF1,DF2] = coefTest(lme,contrasts{iC});
    fprintf('%s: F(%d,%d)=%1.1f, p=%1.2e; ',contrastlabels{iC},DF1,DF2,F,p)
    sigstar(positions{iC},p);
end
writetable(tbl,'SourceData_Fig3d_V1_VarianceExplainedGLM.xlsx')

tbl             = table(Y(idx_AC),X_var(idx_AC),G_mou2(idx_AC),'VariableNames',{'Y','Var','Mouse'}); %Create table for mixed model
lme             = fitlme(tbl,'Y~Var+(1|Mouse)'); %construct linear mixed effects model with fixed effect of temporal window (NO random intercept for mice as cohorts are tested)
stats           = dataset2table(anova(lme,'DFMethod','Satterthwaite')); %Perform ANOVA on model and output as matrix
fprintf('F(%d,%2.0f)=%1.2f, p=%1.2e; \n',stats{2,3},stats{2,4},stats{2,2},stats{2,5})

fprintf('Posthoc comparison:\n')
for iC = 1:length(contrasts)
    [p,F,DF1,DF2] = coefTest(lme,contrasts{iC});
    fprintf('%s: F(%d,%d)=%1.1f, p=%1.2e; ',contrastlabels{iC},DF1,DF2,F,p)
    sigstar(positions{iC}+1,p);
end

for iC = 1:4
    idx = strcmp(X_var,params.Yvarlabels(2+iC));
    tbl             = table(Y(idx),X_area(idx),G_mou2(idx),'VariableNames',{'Y','Area','Mouse'}); %Create table for mixed model
    lme             = fitlme(tbl,'Y~Area+(1|Mouse)'); %construct linear mixed effects model with fixed effect of temporal window (NO random intercept for mice as cohorts are tested)
    stats           = dataset2table(anova(lme,'DFMethod','Satterthwaite')); %Perform ANOVA on model and output as matrix
    fprintf('\n%s: ',params.Yvarlabels{2+iC})
    fprintf('F(%d,%2.0f)=%1.1f, p=%1.2e;',stats{2,3},stats{2,4},stats{2,2},stats{2,5})
    sigstar([1 2] + iC/(params.nShuffleCats+1),stats{2,5});
end
writetable(tbl,'SourceData_Fig3d_AC_VarianceExplainedGLM.xlsx')

set(gca,'XTick',1.5:1:5.5,'XTickLabels',params.areas,'YTick',0:0.02:0.2)
xlim([1 1+params.nAreas])
legend(handles(1:4),{'Trial' 'Visual' 'Audio' 'Motor'},'Location','NorthWest'); legend boxoff;

export_fig(fullfile(params.savedir,sprintf('Box_cvR2_V1_AC_stats')),'-eps','-nocrop')
ylim([-0.01 0.12])
export_fig(fullfile(params.savedir,sprintf('Box_cvR2_V1_AC_axis')),'-eps','-nocrop')

% export_fig(fullfile(params.savedir,sprintf('Bar_cvR2_V1_AC')),'-eps','-nocrop')

%% Computing significance of EV versus shuffling firing rate and prediction:
nShuffles       = 1000;
nNeurons        = length(spikeData.session_ID);

idx_time        = params.xtime>=0 & params.xtime<=0.2e6;

% temp_var_expl_splits = var_expl_splits(expidx,:);
cvR2_shuf         = NaN(nNeurons,nShuffles);
cvR2_cat_cvfit    = NaN(nNeurons,params.nShuffleCats);
cvR2_cat_shuf     = NaN(nNeurons,params.nShuffleCats,nShuffles);
lastsesid           = [];

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
        
        idx_time_all        = repmat(idx_time,1,nTrials); %make an index of all time bins that are included for computation of cvR2

        Y_full              = squeeze(output.y(iNeuron,1:nTotalSpikeBins))'; %get spike rates
        Yh_full             = cvglmnetPredict(output.modelFits(iNeuron,1),X_full,params.lambdastring,'response'); %cv predicted firing rates for full model

        for iC = 1:params.nShuffleCats
            %Predict firing rate when using only one predictor set, compute R2 and store value:
            idx                     = ismember(output.x_label,params.shuffleVars{iC});
            X_temp                  = X_full;
            X_temp(:,~idx)          = 0;
            Yh_cat                  = cvglmnetPredict(output.modelFits(iNeuron,1),X_temp,params.lambdastring,'response');
            
            cvR2_cat_cvfit(iNeuron,iC) = 1 - var(Y_full(idx_time_all) - Yh_cat(idx_time_all)) / var(Y_full(idx_time_all)); %compute R2
            
            for iS = 1:nShuffles
                Y_cat_shuf                = reshape(Yh_cat,params.nTimebins,nTrials); %reshape to time by trial
                Y_cat_shuf                = Y_cat_shuf(:,randperm(nTrials)); %permute these trials, keep NaNs for higher trial numbers the same
                Y_cat_shuf                = reshape(Y_cat_shuf,nTotalSpikeBins,1); %reshape again to vector

                cvR2_cat_shuf(iNeuron,iC,iS) = 1 - var(Y_full(idx_time_all) - Y_cat_shuf(idx_time_all)) / var(Y_full(idx_time_all)); %compute R2
            end
        end
    end
end

%% Determine significance:
alphathr            = 0.001; %Stringent threshold on significance
cvR2_cat_sig        = cvR2_cat > prctile(cvR2_cat_shuf,100-alphathr*100,3);

%% Show explained variance across neurons as bars:
params.labels_cats      = {'Raw' 'Full' 'Trial' 'Vis' 'Aud' 'Motor'};
params.colors_cats      = {[0 0 0] [0.5 0.5 0.5] [0.5 0.2 0.3] [0 0 0.8] [0.8 0 0] [0.1 0.7 0.3]};

idx_exp                 = true(size(spikeData.session_ID)); %all cohorts

idx_area                = strcmp(spikeData.area,'V1');
idx                     = idx_exp & idx_area & any(cvR2_cat_sig(:,2:4),2);

%Sort by aud/vis/motor
temp                    = cvR2_cat_nn(idx,:);
temp                    = temp + abs(randn(size(temp))*1e-3);
% temp = temp+0.01;
sortval                 = (-(temp(:,2) / max(temp(:,2))) + (temp(:,4) / max(temp(:,4)))) ./ (3*temp(:,3) / max(temp(:,3)));
% temp = temp-0.01;
[~,sortidx]             = sort(sortval);
temp                    = temp(sortidx,:); %actual sorting based on derived index

temp_sig                = [];
temp_sig(:,2)           = cvR2_cat_sig(idx,2) & ~cvR2_cat_sig(idx,3) & ~cvR2_cat_sig(idx,4);
temp_sig(:,3)           = ~cvR2_cat_sig(idx,2) & cvR2_cat_sig(idx,3) & ~cvR2_cat_sig(idx,4);
temp_sig(:,4)           = ~cvR2_cat_sig(idx,2) & ~cvR2_cat_sig(idx,3) & cvR2_cat_sig(idx,4);
temp_sig                = temp_sig(sortidx,:); %actual sorting based on derived index

figure; set(gcf,'units','normalized','Position',[0.05 0.59 0.44 0.25],'color','w'); hold all;
for iC = 2:4 %for each of the categories (except trialData)
    subplot(3,1,iC-1); hold all;
    
    h = bar(find(temp_sig(:,iC)==1),temp(temp_sig(:,iC)==1,iC),1.5,'k');
    h.FaceColor = params.colors_cats{iC+2};
    h.EdgeColor = params.colors_cats{iC+2};
    h = bar(find(temp_sig(:,iC)~=1),temp(temp_sig(:,iC)~=1,iC),1.5,'k');
    h.FaceColor = [0.5 0.5 0.5];
%     h.EdgeColor = [1 1 1];
    
    ylim(round(get(gca,'ylim'),1));
    set(gca,'XTick',[],'YTick',get(gca,'ylim'));
    ylabel('cvR2')
    title(params.labels_cats(iC+2))
end

% export_fig(fullfile(params.savedir,sprintf('Bars_cvR2_AV_3cohorts_colored')),'-eps','-nocrop')

%% Fig S6b - explained variance across neurons as scatter based on rank: (stringer et al. 2019 S13)
idx_exp                 = true(size(spikeData.session_ID)); %all cohorts
idx_area                = strcmp(spikeData.area,'V1');
idx                     = idx_exp & idx_area;

I = NaN(sum(idx),4);
for iC = 1:4 %for each of the categories find the ranks
    I(:,iC) = tiedrank(cvR2_cat(idx,iC));
end

temp_sig                = [];
temp_sig(:,2)           = cvR2_cat_sig(idx,2) & ~cvR2_cat_sig(idx,3) & ~cvR2_cat_sig(idx,4);
temp_sig(:,3)           = ~cvR2_cat_sig(idx,2) & cvR2_cat_sig(idx,3) & ~cvR2_cat_sig(idx,4);
temp_sig(:,4)           = ~cvR2_cat_sig(idx,2) & ~cvR2_cat_sig(idx,3) & cvR2_cat_sig(idx,4);

mrkrsize = 6;
figure; set(gcf,'units','normalized','Position',[0.05 0.6 0.37 0.17],'color','w'); hold all;
comparisons = [2 4; 2 3; 3 4];
for iC = 1:3 %for each of the category combinations (except trialData)
    subplot(1,3,iC); hold all;
    
    idx_sig = temp_sig(:,comparisons(iC,1))==0 & temp_sig(:,comparisons(iC,2))==0;
    scatter(I(idx_sig,comparisons(iC,1)),I(idx_sig,comparisons(iC,2)),mrkrsize,'k','filled');
    idx_sig = temp_sig(:,comparisons(iC,1))==1;
    scatter(I(idx_sig,comparisons(iC,1)),I(idx_sig,comparisons(iC,2)),mrkrsize,params.colors_cats{comparisons(iC,1)+2},'filled');
    idx_sig = temp_sig(:,comparisons(iC,2))==1;
    scatter(I(idx_sig,comparisons(iC,1)),I(idx_sig,comparisons(iC,2)),mrkrsize,params.colors_cats{comparisons(iC,2)+2},'filled');

    ylim([0 sum(idx)]);
    xlim([0 sum(idx)]);
    set(gca,'XTick',get(gca,'xlim'),'YTick',get(gca,'ylim'));
    xlabel(sprintf('from %s (rank)',params.labels_cats{comparisons(iC,1)+2}))
    ylabel(sprintf('from %s (rank)',params.labels_cats{comparisons(iC,2)+2}))
end

export_fig(fullfile(params.savedir,sprintf('Scatter_cvR2_cross_AV_StringerS13_colored')),'-eps','-nocrop')

comparisons = [2 4; 2 3; 3 4];
for i = 1:3
    %statistics:
    [r,p] = corr(cvR2_cat(idx,comparisons(i,1)), cvR2_cat(idx,comparisons(i,2)), 'type', 'Spearman','rows' ,'complete');
    fprintf('Spearman rank correlation between %s and %s \n',params.labels_cats{comparisons(i,1)+2},params.labels_cats{comparisons(i,2)+2})
    if p<0.001
        fprintf('r=%1.2f, p=%1.2e; \n',r,p)
    else
        fprintf('r=%1.2f, p=%1.3f; \n',r,p)
    end
end

data        = cvR2_cat(idx,2:4);
data        = [data(:)];
groups      = repmat(params.labels_cats(4:6),sum(idx),1);
groups  	= [groups(:)];
tbl         = table(data,groups,'VariableNames',{'EV','Modality'}); %Create table for mixed model

writetable(tbl,'SourceData_FigS6b_V1_GLMEV_CorrelationAVM.xlsx')

%% Figure with overlap in significant coding in V1:
idx_V1                      = strcmp(spikeData.area,'V1');

params.labels_venn          = {'Vis' 'Aud' 'Mot' 'Vis-Aud'  'Vis-Mot'  'Aud-Mot' 'All'};

params.colors_cats          = {[0 0 0.8] [0.8 0 0] [0.1 0.7 0.3]};

%Make figure:
figure; set(gcf,'units','normalized','Position',[0.08 0.3 0.3 0.42],'color','w')

idx_sign        = cvR2_cat_sig(idx_V1,2:4)';

frac_sign       = NaN(7,1); %store the fraction of neurons significantly coding for this variable

%compute fraction of overlap for each combination:
frac_sign(1)    = sum(idx_sign(1,:) & ~idx_sign(2,:) & ~idx_sign(3,:)) / sum(any(idx_sign));
frac_sign(2)    = sum(~idx_sign(1,:) & idx_sign(2,:) & ~idx_sign(3,:)) / sum(any(idx_sign));
frac_sign(3)    = sum(~idx_sign(1,:) & ~idx_sign(2,:) & idx_sign(3,:)) / sum(any(idx_sign));

frac_sign(4)    = sum(idx_sign(1,:) & idx_sign(2,:) & ~idx_sign(3,:)) / sum(any(idx_sign));
frac_sign(5)    = sum(idx_sign(1,:) & ~idx_sign(2,:) & idx_sign(3,:)) / sum(any(idx_sign));
frac_sign(6)    = sum(~idx_sign(1,:) & idx_sign(2,:) & idx_sign(3,:)) / sum(any(idx_sign));
frac_sign(7)    = sum(idx_sign(1,:) & idx_sign(2,:) & idx_sign(3,:)) / sum(any(idx_sign)); %triple combo

frac_noresp     = sum(~any(idx_sign,1)) / sum(idx_V1);

if round(sum(frac_sign),2)~=1
    error('Does not match to 100%%')
end
frac_sign(frac_sign==0)=0.001;

[H,S] = venn(frac_sign,'FaceColor',params.colors_cats,'FaceAlpha',0.3,'EdgeColor',params.colors_cats);
frac_sign(frac_sign==0.01)=0;

%Now label each zone:
for i = 1:7
    textstring = sprintf('%s %2.0f%%',params.labels_venn{i},frac_sign(i)*100);
    text(S.ZoneCentroid(i,1), S.ZoneCentroid(i,2),textstring)
end
textstring = sprintf('%s %2.0f%%','Non-responsive',frac_noresp*100);
text(-0.7,0.5,textstring)

fprintf('\nPercentage of neurons coding overall:\n\n')
for i = 1:3
    fprintf('%s: %d neurons, %2.1f%%\n',params.labels_venn{i},sum(idx_sign(i,:)),(sum(idx_sign(i,:)) / sum(idx_V1)) *100)
end

%Print results in text
fprintf('\n Total %d neurons:\n\n',sum(idx_V1))
for i = 1:7
    fprintf('%s: %3.0f neurons, %2.1f%%\n',params.labels_venn{i},frac_sign(i)* sum(idx_V1),frac_sign(i)*100)
end
fprintf('%s: %3.0f neurons, %2.1f%%\n','Non-responsive',frac_noresp * sum(idx_V1),frac_noresp*100)

%% Print results in text
z_a = sum(idx_sign(2,:));
z_anm = sum(idx_sign(2,:) & ~idx_sign(3,:));
z_am = sum(idx_sign(2,:) & idx_sign(3,:));
% z_t = z_a +  z_m + z_am;

fprintf('Of the %d auditory neurons, \n',z_a)
fprintf(' %d neurons (%2.1f %%) were exclusively auditory, \n',z_anm, z_anm/z_a*100)
fprintf(' %d neurons (%2.1f %%) were also motor, \n',z_am, z_am/z_a*100)

%% Print results in text
z_a = sum(idx_sign(2,:) & ~idx_sign(3,:));
z_m = sum(~idx_sign(2,:) & idx_sign(3,:));
z_am = sum(idx_sign(2,:) & idx_sign(3,:));
z_t = z_a +  z_m + z_am;

fprintf('Of the %d auditory or motor-responsive neurons, \n',z_t)
fprintf(' %d neurons (%2.1f %%) were exclusively auditory, \n',z_a, z_a/z_t*100)
fprintf(' %d neurons (%2.1f %%) were exclusively motor, \n',z_m, z_m/z_t*100)
fprintf(' %d neurons (%2.1f %%) jointly responsive, \n',z_am, z_am/z_t*100)

%% %%%%%








%% %%%%






%% Compute predicted firing rate response for trial types based on different subpredictors:

fprintf('Computing firing rate response based on subpredictors for neuron        \n');

params.nSplits      = 4;
nNeurons            = length(spikeData.session_ID);
ratemat             = NaN(nNeurons,params.nShuffleCats+2,params.nSplits,params.nTimebins);

lastsesid   = [];

for iNeuron = 1:nNeurons %loop over neurons
    fprintf(repmat('\b', 1, numel([num2str(iNeuron-1) num2str(nNeurons)])+2));
    fprintf('%d/%d\n',iNeuron,nNeurons);
    if ~isempty(output.modelFits(iNeuron,1).lambda)
        
        nTotalSpikeBins     = sum(~isnan(output.y(iNeuron,:)));
        nTrials             = nTotalSpikeBins / params.nTimebins;
        
        temptrialData       = MOL_getTempPerSes(spikeData.session_ID(iNeuron),trialData);

        splits              = {};
        splits{1}           = ismember(temptrialData.trialType,{'X'}) & temptrialData.visualOriChangeNorm==2;
        splits{2}           = ismember(temptrialData.trialType,{'X'}) & temptrialData.visualOriChangeNorm==3;
        
        splits{3}           = ismember(temptrialData.trialType,{'Y'}) & temptrialData.audioFreqChangeNorm==2;
        splits{4}           = ismember(temptrialData.trialType,{'Y'}) & temptrialData.audioFreqChangeNorm==3;
        
        if ~(strcmp(spikeData.session_ID(iNeuron),lastsesid))
            load(fullfile(folderpath,sprintf('X_Full_GLM_%s.mat',spikeData.session_ID{iNeuron})))
            lastsesid           = spikeData.session_ID(iNeuron);
        end
        
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
            idx(1)                  = true; %set trial type to always predict
            X_temp                  = X_full;
            X_temp(:,~idx)          = 0;
            Yh_cat                  = cvglmnetPredict(output.modelFits(iNeuron,1),X_temp,params.lambdastring,'response');
            Yh_cat_r                = reshape(Yh_cat,params.nTimebins,nTrials);

            for iSplit = 1:params.nSplits
                ratemat(iNeuron,iC+2,iSplit,:)     = nanmean(Yh_cat_r(:,splits{iSplit}),2);
            end
        end
    end
end

%% Fig 3f: mean firing rate response to trial types with predictions based on different subsets. 

params.labels_splits    = {'Vthr' 'Vmax' 'Athr' 'Amax'};
params.labels_cats      = {'Raw' 'Full' 'Trial' 'Vis' 'Aud' 'Motor'};
params.colors_cats      = {[0 0 0] [0.5 0.5 0.5] [0.5 0.2 0.3] [0 0 0.8] [0.8 0 0] [0.1 0.7 0.3]};

idx_exp                 = true(size(spikeData.session_ID));

idx_area                = strcmp(spikeData.area,'V1');
idx                     = idx_exp & idx_area;

figure; set(gcf,'units','normalized','Position',[0.25 0.5 0.45 0.28],'color','w'); hold all;

for iSplit = 1:params.nSplits
    subplot(2,params.nSplits,iSplit); hold all;
    for iC = [1 2 4 5 6]
        tmp             = squeeze(nanmean(ratemat(idx,iC,iSplit,:),1));
        plot(params.xtime,tmp,'Color',params.colors_cats{iC},'LineWidth',0.5)
    end
    
    xlim([params.t_pre params.t_post])
    ylim([4.5 9])
    plot([0 0],get(gca,'ylim'),'k:','LineWidth',0.5)
    if iSplit ==4
        set(gca,'XTick',-2e6:0.2e6:2e6,'XTickLabels',(-2e6:0.2e6:2e6)*1e-3, 'YTick',get(gca,'ylim'));
    else
        set(gca,'XTick',-2e6:0.2e6:2e6,'XTickLabels',[],'YTick',[]);
    end
    if iSplit==params.nSplits
        legend(params.labels_cats([1 2 4 5 6]),'FontSize',8);
        legend boxoff
    end
    title(params.labels_splits{iSplit},'FontSize',9)
end

idx_area    = strcmp(spikeData.area,'A1');
idx         = idx_exp & idx_area;

% figure; set(gcf,'units','normalized','Position',[0.05 0.19 0.27 0.28],'color','w'); hold all;
for iSplit = 1:params.nSplits
    subplot(2,params.nSplits,iSplit+params.nSplits); hold all;
    
    for iC = [1 2 4 5 6]
        tmp             = squeeze(nanmean(ratemat(idx,iC,iSplit,:),1));
        plot(params.xtime,tmp,'Color',params.colors_cats{iC},'LineWidth',0.5)
    end
    
    xlim([params.t_pre params.t_post])
    ylim([3 7])
    plot([0 0],get(gca,'ylim'),'k:','LineWidth',0.5)
    if iSplit == 4
        set(gca,'XTick',-2e6:0.2e6:2e6,'XTickLabels',(-2e6:0.2e6:2e6)*1e-3, 'YTick',get(gca,'ylim'));
    else
        set(gca,'XTick',-2e6:0.2e6:2e6,'XTickLabels',[],'YTick',[]);
    end
    
    title(params.labels_splits{iSplit},'FontSize',9)
end

% export_fig(fullfile(params.savedir,sprintf('Rate_overTime_V1_AC_4trialtypes')),'-eps','-nocrop')

%% Fig S6c: same as above but for the three cohorts separately:
% Show mean firing rate response to trial types with predictions based on different subsets

params.labels_splits    = {'Vthr' 'Vmax' 'Athr' 'Amax'};
params.labels_cats      = {'Raw' 'Full' 'Trial' 'Vis' 'Aud' 'Motor'};
params.colors_cats      = {[0 0 0] [0.5 0.5 0.5] [0.5 0.2 0.3] [0 0 0.8] [0.8 0 0] [0.1 0.7 0.3]};

idx_area                = strcmp(spikeData.area,'V1');

figure; set(gcf,'units','normalized','Position',[0.25 0.5 0.45 0.28],'color','w'); hold all;

for iExp = 1:3
    
    idx_exp                 = ismember(spikeData.session_ID,sessionData.session_ID(strcmp(sessionData.Experiment,params.Experiments(iExp))));
    idx                     = idx_exp & idx_area;
    
    for iSplit = [3 4] % Only auditory trials:
        subplot(3,params.nSplits,(iExp-1)*params.nSplits + iSplit); hold all;

        trialoffset = squeeze(nanmean(nanmean(diff(ratemat(idx,[5 3],iSplit,params.xtime<0),[],2),1),4));
        for iC = [1 2 4 5 6]
            tmp             = squeeze(nanmean(ratemat(idx,iC,iSplit,:),1));
            if iC>3 %correct for trial prediction in display:
                tmp = tmp + trialoffset;
            end
            plot(params.xtime,tmp,'Color',params.colors_cats{iC},'LineWidth',0.5)
        end
        xlim([params.t_pre+0.1e6 params.t_post])
        temp = get(gca,'ylim');
        ylim([floor(temp(1)) ceil(temp(2))]);
%             ylim([3 9])
        plot([0 0],get(gca,'ylim'),'k:','LineWidth',0.5)
        set(gca,'XTick',-2e6:0.2e6:2e6,'XTickLabels',(-2e6:0.2e6:2e6)*1e-3, 'YTick',get(gca,'ylim'));
        if iSplit==params.nSplits && iExp==3
            legend(params.labels_cats([1 2 4 5 6]),'FontSize',8);
            legend boxoff
        end
        title(params.labels_splits{iSplit},'FontSize',9)
    end
end

% export_fig(fullfile(params.savedir,sprintf('Rate_overTime_V1_3cohorts_4trialtypes')),'-eps','-nocrop')

%% Fig 3g: For the three cohorts separately
% show auditory and motor response averaged during 0-200 ms

idx_area                = strcmp(spikeData.area,'V1');

response                = squeeze(nanmean(ratemat(:,:,:,params.xtime>0 & params.xtime<=200e3),4)- nanmean(ratemat(:,:,:,params.xtime<0),4));

figure; set(gcf,'units','normalized','Position',[0.25 0.5 0.18 0.29],'color','w'); hold all;

subplot(1,2,1); hold all;
datatotest = squeeze(nanmean(response(idx_area,5,[3 4]),3));

for iExp = 1:3
    idx_exp                 = ismember(spikeData.session_ID,sessionData.session_ID(strcmp(sessionData.Experiment,params.Experiments(iExp))));
    idx                     = idx_area & idx_exp;
    datatoplot              = squeeze(nanmean(response(idx,5,[3 4]),3));
            %         (barplot)
%     h = bar(iExp,nanmean(datatoplot),'FaceColor',params.colors_experiments{iExp});
%     errorbar(iExp,nanmean(datatoplot),nanstd(datatoplot)/sqrt(sum(idx_exp)),'k','LineWidth',1)
    
%     (boxplot)
    h = boxplot(datatoplot, 'plotstyle','compact','positions',iExp,...
    'medianstyle','line','boxstyle','outline','outliersize',0.01,'whisker',1,...
    'colors','k','widths',0.8);
    
end

h = findobj(gca,'tag','Outliers');
delete(h)
h = findobj(gca,'Tag','Box');
for j=1:length(h)
%     patch(get(h(j),'XData'),get(h(j),'YData'),params.colors_experiments{mod(j-1,43)+1},'FaceAlpha',1);
end

ylim([-0.7 2])

tbl             = table(datatotest,X_coh(idx_area),G_mou(idx_area),'VariableNames',{'Y','Var','Mouse'}); %Create table for mixed model
lme             = fitlme(tbl,'Y~Var'); %construct linear mixed effects model with fixed effect of temporal window (NO random intercept for mice as cohorts are tested)
stats           = dataset2table(anova(lme,'DFMethod','Satterthwaite')); %Perform ANOVA on model and output as matrix
fprintf('F(%d,%2.0f)=%1.2f, p=%1.2e; \n',stats{2,3},stats{2,4},stats{2,2},stats{2,5})
sigstar([1 3],stats{2,5})
set(gca,'XTick',1:3,'XTickLabels',params.ExperimentLabels)

writetable(tbl,'SourceData_Fig3g_GLM_AudResponse_V1.xlsx')

% Same but for motor component only:
subplot(1,2,2); hold all;
datatotest = squeeze(nanmean(response(idx_area,6,[3 4]),3));

for iExp = 1:3
    idx_exp                 = ismember(spikeData.session_ID,sessionData.session_ID(strcmp(sessionData.Experiment,params.Experiments(iExp))));
    idx                     = idx_area & idx_exp;
    datatoplot              = squeeze(nanmean(response(idx,6,[3 4]),3));
    %         (barplot)
%     h = bar(iExp,nanmean(datatoplot),'FaceColor',params.colors_experiments{iExp});
%     errorbar(iExp,nanmean(datatoplot),nanstd(datatoplot)/sqrt(sum(idx_exp)),'k','LineWidth',1)
    %     (boxplot)
    h = boxplot(datatoplot, 'plotstyle','compact','positions',iExp,...
    'medianstyle','line','boxstyle','outline','outliersize',0.01,'whisker',1,...
    'colors','k','widths',0.8);
end

h = findobj(gca,'tag','Outliers');
delete(h)
h = findobj(gca,'Tag','Box');
for j=1:length(h)
%     patch(get(h(j),'XData'),get(h(j),'YData'),params.colors_experiments{mod(j-1,43)+1},'FaceAlpha',1);
end

ylim([-0.7 2])

tbl             = table(datatotest,X_coh(idx_area),G_mou(idx_area),'VariableNames',{'Y','Var','Mouse'}); %Create table for mixed model
lme             = fitlme(tbl,'Y~Var'); %construct linear mixed effects model with fixed effect of temporal window (NO random intercept for mice as cohorts are tested)
stats           = dataset2table(anova(lme,'DFMethod','Satterthwaite')); %Perform ANOVA on model and output as matrix
fprintf('F(%d,%2.0f)=%1.2f, p=%1.2e; \n',stats{2,3},stats{2,4},stats{2,2},stats{2,5})

contrasts       = {[0 1 0] [0 0 1] [0 1 -1]};
positions       = {[1 3] [1 2] [2 3]};
contrastlabels  = {'NE vs. MST' 'NE vs. UST' 'UST vs. MST'};

fprintf('Posthoc comparison:\n')
for iC = 1:length(contrasts)
    [p,F,DF1,DF2] = coefTest(lme,contrasts{iC});
    fprintf('%s: F(%d,%d)=%1.1f, p=%1.2e; \n',contrastlabels{iC},DF1,DF2,F,p)
    sigstar(positions{iC},p);
end

ylabel('Evoked response (sp/s)')
set(gca,'XTick',1:3,'XTickLabels',params.ExperimentLabels)

writetable(tbl,'SourceData_Fig3g_GLM_MotorResponse_V1.xlsx')

% export_fig(fullfile(params.savedir,sprintf('Bar_Rate_V1_3Cohorts')),'-eps','-nocrop')
export_fig(fullfile(params.savedir,sprintf('Box_Rate_V1_3Cohorts')),'-eps','-nocrop')


%% Fig S6f: Subsampling analysis:

idx_area                = strcmp(spikeData.area,'V1');
response                = squeeze(nanmean(ratemat(:,:,:,params.xtime>0 & params.xtime<=200e3),4)- nanmean(ratemat(:,:,:,params.xtime<0),4));

iExp = 1;
nSubNeurons = sum(idx_area & ismember(spikeData.session_ID,sessionData.session_ID(strcmp(sessionData.Experiment,params.Experiments(iExp)))));

f = figure; set(gcf,'units','normalized','Position',[0.25 0.5 0.2 0.18],'color','w'); hold all;
handles = [];
subplot(1,2,1); hold all;
datatotest = squeeze(nanmean(response(idx_area,5,[3 4]),3));

for iExp = 1:3
    idx_exp                 = ismember(spikeData.session_ID,sessionData.session_ID(strcmp(sessionData.Experiment,params.Experiments(iExp))));
    idx                     = find(idx_area & idx_exp);
    
    datatoplot  = [];
    for i = 1:1000
        temp = squeeze(nanmean(response(idx(randi(length(idx),nSubNeurons,1)),5,[3 4]),3));
        datatoplot(i,1) = nanmean(temp);
    end
%   barplot:
%     h = bar(iExp,prctile(datatoplot,50),'FaceColor',params.colors_experiments{iExp});
%     errorbar(iExp,prctile(datatoplot,50),prctile(datatoplot,95)-prctile(datatoplot,50),prctile(datatoplot,50)-prctile(datatoplot,5),'k','LineWidth',1)
%   boxplot:
    h = boxplot(datatoplot, 'plotstyle','compact','positions',iExp,...
    'medianstyle','line','boxstyle','outline','outliersize',0.01,'whisker',5,...
    'colors','k','widths',0.8,'datalim',prctile(datatoplot,[5,95]));
end

h = findobj(gca,'tag','Outliers');
delete(h)
h = findobj(gca,'Tag','Box');
for j=1:length(h)
    patch(get(h(j),'XData'),get(h(j),'YData'),params.colors_experiments{mod(j-1,43)+1},'FaceAlpha',1);
end

ylim([0 1.5])

tbl             = table(datatotest,X_coh(idx_area),G_mou(idx_area),'VariableNames',{'Y','Var','Mouse'}); %Create table for mixed model
lme             = fitlme(tbl,'Y~Var'); %construct linear mixed effects model with fixed effect of temporal window (NO random intercept for mice as cohorts are tested)
stats           = dataset2table(anova(lme,'DFMethod','Satterthwaite')); %Perform ANOVA on model and output as matrix
fprintf('F(%d,%2.0f)=%1.2f, p=%1.2e; \n',stats{2,3},stats{2,4},stats{2,2},stats{2,5})
sigstar([1 3],stats{2,5})
set(gca,'XTick',1:3,'XTickLabels',params.ExperimentLabels)

writetable(tbl,'SourceData_FigS6f_GLM_AudResponse_V1.xlsx')

subplot(1,2,2); hold all;
datatotest = squeeze(nanmean(response(idx_area,6,[3 4]),3));

for iExp = 1:3
    idx_exp                 = ismember(spikeData.session_ID,sessionData.session_ID(strcmp(sessionData.Experiment,params.Experiments(iExp))));
    idx                     = find(idx_area & idx_exp);

    datatoplot  = [];
    for i = 1:1000
        temp = squeeze(nanmean(response(idx(randi(length(idx),nSubNeurons,1)),6,[3 4]),3));
        datatoplot(i,1) = nanmean(temp);
    end
%   barplot:
%     h = bar(iExp,prctile(datatoplot,50),'FaceColor',params.colors_experiments{iExp});
%     errorbar(iExp,prctile(datatoplot,50),prctile(datatoplot,95)-prctile(datatoplot,50),prctile(datatoplot,50)-prctile(datatoplot,5),'k','LineWidth',1)
%   boxplot:
    h = boxplot(datatoplot, 'plotstyle','compact','positions',iExp,...
    'medianstyle','line','boxstyle','outline','outliersize',0.01,'whisker',1,...
    'colors','k','widths',0.8,'datalim',prctile(datatoplot,[5,95]));
end

h = findobj(gca,'tag','Outliers');
delete(h)
h = findobj(gca,'Tag','Box');
for j=1:length(h)
    patch(get(h(j),'XData'),get(h(j),'YData'),params.colors_experiments{mod(j-1,43)+1},'FaceAlpha',1);
end

sigstar([1 3],0.009);
sigstar([2 3],0.009);
ylim([0 1.5])

ylabel('Evoked response (sp/s)')
set(gca,'XTick',1:3,'XTickLabels',params.ExperimentLabels)
writetable(tbl,'SourceData_FigS6f_GLM_MotorResponse_V1.xlsx')

% export_fig(fullfile(params.savedir,sprintf('Bar_Rate_V1_3Cohorts_subsampled')),'-eps','-nocrop')
export_fig(fullfile(params.savedir,sprintf('Box_Rate_V1_3Cohorts_subsampled')),'-eps','-nocrop')

%% %





%% As a function of depth:


%% Depth corrections based on histology:
tempcomb = {'20212019070104' 150; '20232019070170' -100; '20282019121322' -100; '20282019121626' -100; '20292019121240' -50; '20292019121343' -50};

for iSes = 1:size(tempcomb,1)
    spikeData.ChannelY(strcmp(spikeData.session_ID,tempcomb{iSes,1})) = spikeData.ChannelY(strcmp(spikeData.session_ID,tempcomb{iSes,1})) + tempcomb{iSes,2};
end

%%
params.depthcorrection          = -50;

spikeData.ChannelY              = spikeData.ChannelY + params.depthcorrection;

%% Fig 4: Make SUA heatmap over laminar depth for each trial type and for each model split:

%Histogram with binning on depth:
params.chdepthmin       = -12.5;
params.chdepthmax       = 1012.5;
params.resolution       = 30;

binedges                = params.chdepthmin:params.resolution:params.chdepthmax;
nBins                   = length(binedges)-1;

params.finalYaxis       = binedges(1:end-1)+params.resolution/2;

nModelSplits            = size(ratemat,2);
SUAdepthmap_all         = zeros(nModelSplits,params.nSplits,nBins,params.nTimebins);

idx_area                = strcmp(spikeData.area,'V1');

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
spat_filter     = fspecial('average',[11 4]);

for iMod = 1:nModelSplits
    for iSplit = 1:params.nSplits
        temp = squeeze(SUAdepthmap_all(iMod,iSplit,:,:));
        SUAdepthmap_all(iMod,iSplit,:,:) = conv2(temp,spat_filter,'same');
    end
end

%
plotModels      = [1 4 6 1 5 6];
plotTrialTypes  = {[1 2] [1 2] [1 2] [3 4] [3 4] [3 4]};

figure; set(gcf,'units','normalized','Position',[0.05 0.5 0.4 0.3],'color','w');
timeticks = [-0.2e6 0 0.2e6 0.4e6 0.6e6];
binticks = [0 250 500 750 1000];

cmapdata = parula(128);

for iMod = 1:length(plotModels)
    subplot(2,3,iMod);        hold all;
    
    title(params.labels_cats(plotModels(iMod)),'FontSize',10);
    
    SUAdepthmap = squeeze(nanmean(SUAdepthmap_all(plotModels(iMod),plotTrialTypes{iMod},:,:),2));
    
    imagesc(params.xtime,params.finalYaxis,SUAdepthmap); 
    plot([0 0],[-2000 2000],':','Color',[1 1 1],'LineWidth',1)
    temp = caxis;
    caxis([-0.3 temp(2)])
    colormap(cmapdata)

    %Figure make up:
    if iMod==1 || iMod==4
        ylabel('Cortical depth (\mum)','FontSize', 9)
    end
    if iMod==5
        xlabel('Time (ms)','FontSize',9)
    end
    set(gca,'YDir','reverse') %correct for imagesc plotting inverse
    set(gca,'XTick',timeticks,'XTickLabel',timeticks*1e-3)
    set(gca,'YTick',binticks,'YTickLabel',binticks,'fontsize',7,'tickdir','out')
    xlim([find(params.xtime>-0.2e6,1) find(params.xtime<0.6e6,1,'last')])
    xlim([-0.2e6 0.6e6]);
    ylim([0 1000])
    
    if iMod==3
        temppos = get(gca,'position');
        colorbar('Location','eastoutside');
        set(gca,'Position',temppos)
    end
end

export_fig(fullfile(params.savedir,sprintf('LaminarHeatmap_ModelSplits_AV_3cohorts')),'-eps','-nocrop')


%% Fig S6d,e Control figure, volume vs cvR2 
params.colors_cats      = {[0.5 0.2 0.3] [0 0 0.8] [0.8 0 0] [0.1 0.7 0.3]};

spikeData.dblSoundLevel = NaN(nNeurons,1);
for iSes = 1:nSessions
    spikeData.dblSoundLevel(strcmp(spikeData.session_ID,sessionData.session_ID(iSes))) = sessionData.dblSoundLevel(iSes);
end

figure; set(gcf,'units','normalized','Position',[0.2 0.43 0.22 0.3],'color','w'); hold all;
subplot(2,1,1); hold all;
h = scatter(spikeData.dblSoundLevel+randn(size(spikeData.dblSoundLevel))*1e-3,cvR2_cat(:,3),15,params.colors_cats{3},'.');
xlim([0.5 0.85])
ylim([0 0.4])
ylabel('Aud. EV')
[r,p]               = corr(spikeData.dblSoundLevel,cvR2_cat(:,3),'rows','complete','type', 'Spearman');
fprintf('Spearman rank correlation between sound volume \n and auditory-related EV, r=%1.3f, p=%1.3f; \n',r,p)
subplot(2,1,2); hold all;
h = scatter(spikeData.dblSoundLevel+randn(size(spikeData.dblSoundLevel))*1e-3,cvR2_cat(:,4),15,params.colors_cats{4},'.');
xlim([0.5 0.85])
ylim([0 0.4])
ylabel('Motor EV')
xlabel('Session sound volume (dBA)')

[r,p]               = corr(spikeData.dblSoundLevel,cvR2_cat(:,4),'rows','complete','type', 'Spearman');
fprintf('Spearman rank correlation between sound volume \n and motor-related EV, r=%1.3f, p=%1.3f; \n',r,p)
export_fig(fullfile(params.savedir,sprintf('EV_vs_soundvolume_scatter_Aud_Motor')),'-eps','-nocrop')

tbl             = table(cvR2_cat(:,3),cvR2_cat(:,4),spikeData.dblSoundLevel,...
    'VariableNames',{'Aud EV','Motor EV','Volume'}); %Create table for mixed model
writetable(tbl,'SourceData_FigS6de_GLM_VolumeControl.xlsx')

