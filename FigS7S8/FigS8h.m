%% Oude Lohuis et al. 2023 Nat Neurosci
% Triple dissociation of auditory, visual, and motor processing in primary visual cortex
% MOL (C) 2023

% This script computes the onset latency of neuron during auditory or visual trials 
% an audiovisual change detection task for three different task contingencies 
% as reported in Oude Lohuis et al. 2023

%%
startover

%% Parameter settings:

params.Experiments          = {'ChangeDetectionConflictDecor' 'VisOnlyTwolevels' 'ChangeDetectionConflict'}; %Which versions of the task to load data from
params.ExperimentLabels     = {'NE' 'UST' 'MST'}; %Which versions of the task to load data from
params.AlignOn              = 'stimChange';         %On which timestamp to align trials to
params.nExperiments         = length(params.Experiments);
params                      = MOL_getColors_CHDET(params);

params                      = params_histresponse_auV1(params);

params.areas                = {'V1'};
params.nAreas               = length(params.areas);

params.savedir              = 'E:\Documents\PhD\Figures\Project CHDET\Results - auV1\5LaminarCircuit';

%% Get data:
[Data] = MOL_GetData('E:','CHDET',params.Experiments,{},[],{'sessionData' 'trialData_newtrials' 'spikeData'});
sessionData         = Data.sessionData;
trialData           = Data.trialData;
spikeData           = Data.spikeData;

%% Remove last 20 trials:
trialData           = MOL_RemoveLastnTrials(trialData,20);

%% Filter out neurons based on quality:
spikeData           = MOL_filterNeurons(sessionData,trialData,spikeData);

%% Remove sessions that are passive/active from 2044 and 2045s: 
sesids                      = sessionData.session_ID(~ismember(sessionData.mousename,{'2044' '2045'}));
fprintf('Removed %d/%d sessions from 2044 and 2045 with active and passive epochs\n',length(sessionData.session_ID)-length(sesids),length(sessionData.session_ID));
[sessionData, trialData, spikeData]        = MOL_getTempPerSes(sesids,sessionData,trialData,spikeData);

% %% Remove sessions that are passive: 
% sesids                      = sessionData.session_ID(strcmp(sessionData.State,'Behaving'));
% fprintf('Removed %d/%d sessions with passive stimulation\n',length(sessionData.session_ID)-length(sesids),length(sessionData.session_ID));
% [sessionData, trialData, spikeData]        = MOL_getTempPerSes(sesids,sessionData,trialData,spikeData);

%% Filter out sessions with muscimol:
sesids              = sessionData.session_ID(cellfun(@isempty,sessionData.MuscimolArea));
fprintf('Removed %d/%d sessions  with muscimol manipulations\n',length(sessionData.session_ID)-length(sesids),length(sessionData.session_ID));
[sessionData, trialData, spikeData]        = MOL_getTempPerSes(sesids,sessionData,trialData,spikeData);

%% Filter out trials with photostimulation in V1:
sesids              = sessionData.session_ID(sessionData.UseOpto & strcmp(sessionData.PhotostimArea,'V1'));
idx                 = trialData.hasphotostim==1 & ismember(trialData.session_ID,sesids);
trialFields         = fieldnames(trialData);
fprintf('Removed %d/%d trials  with photostim in V1\n',sum(idx),length(trialData.session_ID));
for iF = 1:length(trialFields)
    trialData.(trialFields{iF}) = trialData.(trialFields{iF})(~idx);
end

%% Filter neurons on area:
spikeFields     = fieldnames(spikeData);
idx             = ismember(spikeData.area,params.areas);
fprintf('Filtered %d/%d neurons based on area\n',sum(idx),length(spikeData.session_ID));
for iF = 1:length(spikeFields)
    spikeData.(spikeFields{iF}) = spikeData.(spikeFields{iF})(idx);
end

%% Get only sessions with neurons recorded:
sesids              = unique(spikeData.session_ID);
fprintf('Removed %d/%d sessions without neurons in target area\n',length(sessionData.session_ID)-length(sesids),length(sessionData.session_ID));
[sessionData,trialData,spikeData] = MOL_getTempPerSes(sesids,sessionData,trialData,spikeData);

%% Report dataset:
fprintf('Dataset: %d sessions, %d trials, %d neurons\n',length(sessionData.session_ID),length(trialData.session_ID),length(spikeData.session_ID));

%%
load('E:\Matlab\Dataset1.mat')
params.savedir              = 'E:\Documents\PhD\Figures\Project CHDET\Results - auV1\5LaminarCircuit';

% %Get only sessions with correct postchange freqs and oris:
% sesids              = unique(trialData.session_ID(~isnan(trialData.visualOriPostChangeNorm) & ~isnan(trialData.audioFreqPostChangeNorm)));
% [sessionData,trialData,spikeData] = MOL_getTempPerSes(sesids,sessionData,trialData,spikeData);

fprintf('Dataset: %d sessions, %d trials, %d neurons\n',length(sessionData.session_ID),length(trialData.session_ID),length(spikeData.session_ID));








%% Compute latencies:

nNeurons                = length(spikeData.ts);
lastsesid               = []; %memory var to keep track of whether neuron comes from different session and design matrix needs to be reconstructed
fprintf('Computing ZETA for neuron        \n');

params.nSplits          = 4;
matLat                  = NaN(nNeurons,params.nSplits);
matP                    = NaN(nNeurons,params.nSplits);
matZETA                 = NaN(nNeurons,params.nSplits);

params.nTimebins_Zeta   = 100; %interpolate to same time axis:

matRate                 = NaN(nNeurons,params.nSplits,params.nTimebins_Zeta);

params.t_pre            = -0.2e6;
params.minTrialCond     = 10;

%Loop over all neurons and compute latency to spike:
for iNeuron = 1:nNeurons %Loop over all neurons:
    fprintf(repmat('\b', 1, numel([num2str(iNeuron-1) num2str(nNeurons)])+2));
    fprintf('%d/%d\n',iNeuron,nNeurons);

    if ~strcmp(lastsesid,spikeData.session_ID(iNeuron)) %construct new predictor matrix if neuron comes from a new session:
        %Get the relevant data for each session individually:
        temptrialData        = MOL_getTempPerSes(spikeData.session_ID(iNeuron),trialData);
        lastsesid            = spikeData.session_ID(iNeuron); %save this session_ID
    end

    %Compute histogram:
    events_ts               = temptrialData.(params.AlignOn);

    splits                  = {};
    splits{1}               = ismember(temptrialData.trialType,{'X'}) & temptrialData.visualOriChangeNorm==3 & ismember(temptrialData.visualOriPostChangeNorm,[1 2]);
    splits{2}               = ismember(temptrialData.trialType,{'X'}) & temptrialData.visualOriChangeNorm==3 & ismember(temptrialData.visualOriPostChangeNorm,[3 4]);
    splits{3}               = ismember(temptrialData.trialType,{'Y'}) & temptrialData.audioFreqChangeNorm==3 & ismember(temptrialData.audioFreqPostChangeNorm,[1 2]);
    splits{4}               = ismember(temptrialData.trialType,{'Y'}) & temptrialData.audioFreqChangeNorm==3 & ismember(temptrialData.audioFreqPostChangeNorm,[3 4]);

    respwin                     = [0 -params.t_pre+0.2e6]*1e-6;

    for iSplit = 1:params.nSplits
        if sum(splits{iSplit})>=params.minTrialCond

            [~,templat,tempz,tempr]     = getZeta(spikeData.ts{iNeuron}*1e-6,(events_ts(splits{iSplit})+params.t_pre)*1e-6,1,[],0,4,respwin,[]);
            if ~all(isnan(templat))
                matLat(iNeuron,iSplit)      = templat(4) + params.t_pre*1e-6;
                matP(iNeuron,iSplit)        = tempz.dblP;
                matZETA(iNeuron,iSplit)     = tempz.dblZETA;

                %time resolution can get too high for double, so remove unique from
                [~,uidx] = unique(tempr.vecT); tempr.vecT = tempr.vecT(uidx); tempr.vecRate = tempr.vecRate(uidx);

                matRate(iNeuron,iSplit,:)   = interp1(tempr.vecT,tempr.vecRate,linspace(0,1,params.nTimebins_Zeta));
            end
        end
    end
end

%% init fields:

spikeData.sign_visresponse      = any(matP(:,[1 2])<0.05,2);
spikeData.sign_auresponse       = any(matP(:,[3 4])<0.05,2);

spikeData.z_visresponse         = max(matZETA(:,[1 2]),[],2);
spikeData.z_auresponse          = max(matZETA(:,[3 4]),[],2);

%% get latencies:
vecLat_vis = NaN(nNeurons,1);
vecLat_aud = NaN(nNeurons,1);

for iNeuron = 1:nNeurons
    [maxZ,idx]            = max(matZETA(iNeuron,[1 2]));
    if maxZ>2
        vecLat_vis(iNeuron,1)      = matLat(iNeuron,idx);
    end
    [maxZ,idx]            = max(matZETA(iNeuron,[3 4]));
    if maxZ>2
        vecLat_aud(iNeuron,1)      = matLat(iNeuron,idx+2);
    end
    
end

vecLat_vis                  = vecLat_vis*1e6;
vecLat_aud                  = vecLat_aud*1e6;

vecLat_vis(vecLat_vis<0)    = NaN;
vecLat_aud(vecLat_aud<0)    = NaN;

%%
figure; hold all; set(gcf,'units','normalized','Position',[0.45 0.5 0.24 0.23],'color','w');

binres      = 10e3;
binedges    = -30e3:binres:200e3;
bintimes    = binedges(2:end);

tempylim = [0 40];
tempxlim = [0 160e3];

% iExp = 2;
% idx_exp = ismember(spikeData.session_ID,sessionData.session_ID(strcmp(sessionData.Experiment,params.Experiments(iExp))));
idx_exp = true(size(spikeData.session_ID));
idx_V1  = strcmp(spikeData.area,'V1');

hcounts = histcounts(vecLat_vis(spikeData.sign_visresponse & idx_V1 & idx_exp),binedges,'Normalization','count');
plot(bintimes,hcounts,'b','LineWidth',1);
plot(nanmedian(vecLat_vis(spikeData.sign_visresponse & idx_V1 & idx_exp)),tempylim(2)-5,'bv','MarkerSize',12,'LineWidth',1);

hcounts = histcounts(vecLat_aud(spikeData.sign_auresponse & idx_V1 & idx_exp),binedges,'Normalization','count');
plot(bintimes,hcounts,'r','LineWidth',1);
plot(nanmedian(vecLat_aud(spikeData.sign_auresponse & idx_V1 & idx_exp)),tempylim(2)-5,'rv','MarkerSize',12,'LineWidth',1);

xlim(tempxlim)
ylim(tempylim)
set(gca,'XTick',0:20e3:160e3,'XTickLabel',(0:20:160),'YTick',tempylim(2))
ylabel('Neuron count')
xlabel('Onset latency (ms)')
title('V1')

%Statistics:
G_mou       = cell(nNeurons,1);
uMice       = unique(sessionData.mousename);
for iMouse = 1:length(uMice)
    G_mou(ismember(spikeData.session_ID,sessionData.session_ID(strcmp(sessionData.mousename,uMice{iMouse})))) = uMice(iMouse);
end

idx_v       = spikeData.sign_visresponse & idx_V1;
idx_a       = spikeData.sign_auresponse & idx_V1;
Y           = [vecLat_vis(idx_v); vecLat_aud(idx_a)];
X_mod       = [ones(sum(idx_v),1); ones(sum(idx_a),1)*2];
X_mou       = [G_mou(idx_v); G_mou(idx_a)]; 
tbl         = table(Y,X_mod,X_mou,'VariableNames',{'Latency','Modality','Mouse'}); %Create table for mixed model

lme             = fitlme(tbl,'Latency~Modality+(1|Mouse)'); %construct linear mixed effects model with fixed effect of temporal window and random intercept for different mice
stats           = dataset2table(anova(lme,'DFMethod','Satterthwaite')); %Perform ANOVA on model and output as matrix
p = stats{2,5};
if p<0.05
    sigstar([nanmedian(vecLat_vis(spikeData.sign_visresponse & idx_V1 & idx_exp)) nanmedian(vecLat_aud(spikeData.sign_auresponse & idx_V1 & idx_exp))],p)
end

%statistics: latencies
fprintf('Spiking onset was significantly earlier for auditory versus visual stimuli:\n')
fprintf('%3.1f ms (%3.1f - %3.1f ms) vs %3.1f ms (%3.1f - %3.1f ms); \n',prctile(vecLat_vis(idx_v),[50 25 75])*1e-3,prctile(vecLat_aud(idx_a),[50 25 75])*1e-3)
fprintf('F(%d,%2.0f)=%1.2f, p=%1.3e; \n',stats{2,3},stats{2,4},stats{2,2},stats{2,5})

filename = sprintf('Hist_Lat_V1_AV.eps');
export_fig(fullfile(params.savedir,filename),gcf);

%% Laminar profile:
%depthcorrections:
tempcomb = {'20212019070104' 150; '20232019070170' -100; '20282019121322' -100; '20282019121626' -100; '20292019121240' -50; '20292019121343' -50};
for iSes = 1:size(tempcomb,1)
    spikeData.ChannelY(strcmp(spikeData.session_ID,tempcomb{iSes,1})) = spikeData.ChannelY(strcmp(spikeData.session_ID,tempcomb{iSes,1})) + tempcomb{iSes,2};
end
params.depthcorrection = -50;
spikeData.ChannelY = spikeData.ChannelY + params.depthcorrection;

%% Figure for layer figure:
params.markersize       = 25;

figure; hold all; set(gcf,'units','normalized','Position',[0.55 0.45 0.27 0.38],'color','w');
idx_area    = strcmp(spikeData.area,'V1');
idx         = spikeData.sign_visresponse & idx_area & spikeData.celltype==1;
scatter(vecLat_vis(idx),-spikeData.ChannelY(idx),(spikeData.z_visresponse(idx)-1)*10,'b^','filled')
idx         = spikeData.sign_visresponse & idx_area & spikeData.celltype==2;
scatter(vecLat_vis(idx),-spikeData.ChannelY(idx),(spikeData.z_visresponse(idx)-1)*10,'bo','filled')

idx         = spikeData.sign_auresponse & idx_area & spikeData.celltype==1;
scatter(vecLat_aud(idx),-spikeData.ChannelY(idx),(spikeData.z_auresponse(idx)-1)*10,'r^','filled')
idx         = spikeData.sign_auresponse & idx_area & spikeData.celltype==2;
scatter(vecLat_aud(idx),-spikeData.ChannelY(idx),(spikeData.z_auresponse(idx)-1)*10,'ro','filled')

%Plot mean of bins:
binres                  = 75;
binedges                = 0:binres:1000;
nBins                   = length(binedges)-1;
depthYaxis              = binedges(2:end)-binres/2; %correct edges

latY_visual             = NaN(nBins,1); %Init variable
latY_audio              = NaN(nBins,1); %Init variable

for iBin = 1:nBins
    idx = spikeData.sign_visresponse & idx_area & spikeData.ChannelY>=binedges(iBin) & spikeData.ChannelY<binedges(iBin+1);
    if sum(idx)>2
        latY_visual(iBin) = nanmedian(vecLat_vis(idx),1);
    end
    idx = spikeData.sign_auresponse & idx_area & spikeData.ChannelY>=binedges(iBin) & spikeData.ChannelY<binedges(iBin+1);
    if sum(idx)>2
        latY_audio(iBin) = nanmedian(vecLat_aud(idx & vecLat_aud<100e3),1); %#ok<*NANMEDIAN>
%         latY_audio(iBin) = nanmedian(vecLat_aud(idx),1); %#ok<*NANMEDIAN>
    end
end

plot(latY_visual,-depthYaxis,'-b','LineWidth',2)
plot(latY_audio,-depthYaxis,'-r','LineWidth',2)

[r,p] = corr(latY_visual,depthYaxis','rows','complete');
fprintf('Cortical depth was significantly correlated to spiking onset latency during visual trials (r=%1.2f, p=%1.3f)\n',r,p)
text(110e3,-100,sprintf('*\n r=%1.2f\n p=%1.3f\n',r,p),'Color',[0 0 1])

[r,p] = corr(latY_audio,depthYaxis','rows','complete');
fprintf('but not auditory trials (r=%1.2f, p=%1.3f)\n',r,p)
text(40e3,-100,sprintf('n.s.\n r=%1.2f\n p=%1.3f\n',r,p),'Color',[1 0 0])

xlim([0 160e3])
ylim([-1050 0])
set(gca,'XTick',0:20e3:160e3,'XTickLabel',(0:20:160),'YTick',-1000:250:0)
ylabel('Cortical depth (\mum)')
xlabel('Onset latency')
MOL_prepfigAI

filename = sprintf('LaminarHist_PvPyr_V1_AV.eps');
export_fig(fullfile(params.savedir,filename),gcf);



%% END OF SCRIPT