%% Oude Lohuis et al. 2023 Nat Neurosci
% Triple dissociation of auditory, visual, and motor processing in primary visual cortex
% MOL (C) 2023

% Script analyzes the GLM results of fitting increasing numbers of video PCs to explain 
% firing rate in V1

startover

%% Load dataset:
folderpath                  = fullfile('E:','Matlab','oudelohuis-et-al-2023-natneurosci','FigS6','GLMfits_motor');
fileList                    = dir(fullfile(folderpath,'*.mat'));
fileList                    = {fileList(:).name};
fileList                    = fileList(~contains(fileList,'X'));
nFiles                      = length(fileList);

%Load the first session, then concatenate the rest
loadstruct          = load(fullfile(folderpath,fileList{1}));

output.x_sesid      = cell(nFiles,1);
output.y            = NaN(2000,size(loadstruct.output.y,2));
output.modelFits    = loadstruct.output.modelFits;
output.x_label      = loadstruct.output.x_label;
nTotalneurons       = 0;

sessionData         = struct();
trialData           = struct();
spikeData           = struct();

cvR2                = [];
cvR2_motor          = [];

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
    
    cvR2            = [cvR2; loadstruct.cvR2]; %#ok<*AGROW>
    cvR2_motor      = [cvR2_motor; loadstruct.cvR2_motor];
    
    nTotalneurons       = nTotalneurons+nNeurons; 
    
end

sessionData.Experiment      = strrep(sessionData.Experiment,num2str(2),'');
params.savedir              = 'E:\Documents\PhD\Figures\Project CHDET\Results - auV1\18GLM\';
cvR2_motor(cvR2_motor<0)    = NaN;%Correct for negative R2

%% Fig S6a: 
% Show cvR2 for consecutive video PCs:
params.nSVDs = 128;

params.areas            = {'V1'};% 'PPC' 'CG1'};
params.nAreas           = length(params.areas); 

figure; set(gcf,'units','normalized','Position',[0.2 0.43 0.22 0.3],'color','w'); hold all;

params.colors_splits = {[0.4 0.4 0.4] [0 0 0.8] [0.8 0 0] [0.1 0.7 0.3]};
handles = [];
idx_V1     = strcmp(spikeData.area,'V1');
h = shadedErrorBar(1:params.nSVDs,nanmean(cvR2_motor(idx_V1,:),1),nanstd(cvR2_motor(idx_V1,:),[],1) / sqrt(nNeurons),{'-','Color',[0 0 0],'LineWidth',2},0);
delete(h.edge(:));

ylabel('Explained Variance')
xlabel('Video PC #')
xlim([0.8 params.nSVDs])
ylim([0 0.035])
set(gca,'XTick',[1 4 16 64 128],'YTick',[0 0.005 0.01 0.015 0.02])
set(gca,'XTick',[1 5 10 20 50 128],'YTick',[0 0.01 0.02 0.03])
set(gca,'xscale','log')

export_fig(fullfile(params.savedir,sprintf('cvR2_motordim')),'-eps','-nocrop')
