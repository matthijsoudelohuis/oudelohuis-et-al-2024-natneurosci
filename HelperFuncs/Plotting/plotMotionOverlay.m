function plotMotionOverlay(params)

params.rawdatadir       = fullfile(params.rootDir,'Data','CHDET','RawData',params.animal,params.session,params.experiment);

videoFile               = fullfile(params.rawdatadir,[params.videoName '.mpg']); % 
csvFile                 = fullfile(params.rawdatadir,[params.videoName 'DLC_resnet50_PupilTrackingMay29shuffle1_1030000.csv']); %
alignFile               = fullfile(params.rawdatadir,[params.videoName '_videoalign.mat']); 
procFile                = fullfile(params.rawdatadir,[params.videoName '_proc.mat']); %
[Data]                  = MOL_GetData('E:','CHDET',{'ChangeDetectionConflict' 'ChangeDetectionConflictDecor' 'VisOnlyTwolevels'},{},params.session,{'sessionData' 'trialData_newtrials' 'spikeData' 'videoData'});

%% Load processes video data
if ~exist(procFile,'file')
    warning('No processed video data file found ')
end
if ~exist(alignFile,'file')
    warning('No mat file with alignmentdata')
end
load(alignFile)
load(procFile)

%% load sessionData, trialData, videoData
sessionData     = Data.sessionData;
trialData       = Data.trialData;
spikeData       = Data.spikeData;
videoData       = Data.videoData;

%% Remove last 5 trials:
trialData = MOL_RemoveLastnTrials(trialData,5);

%% Initialize video reading:
cap                     = py.cv2.VideoCapture(videoFile);

%Read properties:
fps                     = cap.get(int8(5)); %#https://docs.opencv.org/2.4/modules/highgui/doc/reading_and_writing_images_and_video.html#videocapture-get
nframes                 = cap.get(int8(7));
widthPixels             = cap.get(int8(3));
heightPixels            = cap.get(int8(4));

%% try to read some frames until it works, then reset: 
while cap.isOpened
    ret     = cap.read;
    mt      = cell(ret(1));
    if mt{1}
        iFr=1;
        break;
    end
end
cap.set(int8(1),0);

%% Figure 
figure; set(gcf,'units','normalized','Position',[0.45 0.4 0.24 0.46],'color','w')
temp = reshape(proc.avgframe,size(proc.wpix{1},1),size(proc.wpix{1},2));
imagesc(flipud(temp'))
colormap(gray)
ylabel(''); xlabel(''); set(gca,'XTick',[],'YTick',[]);

%% 
figure; set(gcf,'units','normalized','Position',[0.45 0.4 0.24 0.46],'color','w')

temp = reshape(proc.avgmotion,size(proc.wpix{1},1),size(proc.wpix{1},2));

imagesc(flipud(temp'))

h       =  coolwarm(256)  + repmat(0.1-abs(linspace(-0.1,.1,256))',1,3);
h(h>1)  = 1;

% h = getPyPlot_cMap('seismic');

colormap(h);

% lim = max(abs(min(temp(:))),max(temp(:)));
% caxis([-lim lim]);
    
ylabel(''); xlabel(''); set(gca,'XTick',[],'YTick',[]);

%% 
% dims = round(proc.ROI{1}{1}([3 4])) + 1;
dims = round(proc.ROI{1}{1}([3 4]));
nSVDs = 10;

h = getPyPlot_cMap('seismic');

figure; set(gcf,'units','normalized','Position',[0.05 0.4 0.58 0.28],'color','w')
for iSVD = 1:nSVDs
    subplot(1,nSVDs,iSVD)
    temp = reshape(proc.uMotMask{1}(:,iSVD),dims(2),dims(1));
    imagesc(flipud(temp'))
    
%     h       =coolwarm(256) + repmat(0.1-abs(linspace(-0.1,.1,256))',1,3);
%     h(h>1)  = 1;
    
%     h       = coolwarm(256) + repmat(0.1-abs(linspace(-0.1,.1,256))',1,3);
%     h(h>1)  = 1;
    lim = max(abs(min(temp(:))),max(temp(:)));
    caxis([-lim lim]);
    
    colormap(h);
    ylabel(''); xlabel(''); set(gca,'XTick',[],'YTick',[]);
    title(sprintf('PCA %d',iSVD));
end


%% Shift all motion SVD values by frames given that motion is computed as diff between frames
nSessions = length(sessionData.session_ID);
fprintf('Shifting frames for session      \n');
for iSes = 1:nSessions
    fprintf(repmat('\b', 1, numel([num2str(iSes-1) num2str(nSessions)])+1));
    fprintf('%d/%d',iSes,nSessions);
%     videoData.motSVD{iSes} = [nan(1,500); videoData.motSVD{iSes}(1:end-1,:)];
    videoData.motSVD{iSes} = circshift(videoData.motSVD{iSes},params.shiftNframes,1);
end
fprintf('\n')

%% 

params.AlignOn              = 'stimChange';      %On which timestamp to align as t=0

params.videofield           = 'motSVD';

params.t_pre                = -0.5e6;
params.t_post               = 2e6;

params.fs                   = 25; %Hz

params.nSVDs                = 500;

[hist_mat,hist_mat_tot,hist_mat_z,params] = MOL_PSTH_video(sessionData,trialData,videoData,params);

% hist_mat_z_svd              = hist_mat;
% 
% for iSVD = 1:nSVDs
%     temp                        = hist_mat_z_svd(params.xtime_video<0,iSVD,:);
%     hist_mat_z_svd(:,iSVD,:)    = (hist_mat_z_svd(:,iSVD,:) - repmat(nanmean(temp,1),params.nTimebins_video,1,1)) / nanstd(temp(:));
% end

% %%
% 
% hist_mat            = NaN(params.nTimebins_video,nSVDs,nTrials);
% 
% % temp                = zscore(videoData.(params.videofield){ses_idx},[],2);
% 
% fprintf('Computing svd response for trial        \n');
% for iTrial = 1:nTrials
%     fprintf(repmat('\b', 1, numel([num2str(iTrial-1) num2str(nTrials)])+1));
%     fprintf('%d/%d',iTrial,nTrials);
%     ses_idx                             = strcmp(sessionData.session_ID,trialData.session_ID(iTrial));
%     idx                                 = videoData.ts{ses_idx}>trialData.(params.AlignOn)(iTrial)+params.t_pre & videoData.ts{ses_idx}<trialData.(params.AlignOn)(iTrial)+params.t_post;
% %     hist_mat(1:sum(idx),:,iTrial)     	= temp(idx,1:nSVDs);
%     hist_mat(1:sum(idx),:,iTrial)     	= videoData.(params.videofield){ses_idx}(idx,1:nSVDs);
% end
% fprintf('\n')
% 
% %% 
% % hist_mat_z  = NaN(params.nTimebins_video,nTrials);
% nSessions                       = length(sessionData.session_ID);
% 
% hist_mat_tot                    = squeeze(nansum(abs(hist_mat(:,1:nSVDs,:)),2));
% hist_mat_z                      = hist_mat;
% 
% for iSVD = 1:nSVDs
%     temp                        = hist_mat_z(params.xtime_video<0,iSVD,:);
%     hist_mat_z(:,iSVD,:)        = (hist_mat_z(:,iSVD,:) - repmat(nanmean(temp,1),params.nTimebins_video,1,1)) / nanstd(temp(:));
% end

%% 

% nSplits = 1;

% splits          = {};
% splits{1}       = strcmp(trialData.trialType,'X') & trialData.visualOriChangeNorm==2;
% splits{2}       = strcmp(trialData.trialType,'X') & trialData.visualOriChangeNorm==2;
% splits{3}       = strcmp(trialData.trialType,'Y') & trialData.audioFreqChangeNorm==3;
% splits{4}       = strcmp(trialData.trialType,'Y') & trialData.audioFreqChangeNorm==3;

splits          = {};
splits{1}       = strcmp(trialData.trialType,'X') & trialData.visualOriChangeNorm==2 & ismember(trialData.visualOriPostChangeNorm,[1 2]);
splits{2}       = strcmp(trialData.trialType,'X') & trialData.visualOriChangeNorm==2 & ismember(trialData.visualOriPostChangeNorm,[3 4]);
splits{3}       = strcmp(trialData.trialType,'X') & trialData.visualOriChangeNorm==3 & ismember(trialData.visualOriPostChangeNorm,[1 2]);
splits{4}       = strcmp(trialData.trialType,'X') & trialData.visualOriChangeNorm==3 & ismember(trialData.visualOriPostChangeNorm,[3 4]);
splits{5}       = strcmp(trialData.trialType,'Y') & trialData.audioFreqChangeNorm==2 & ismember(trialData.audioFreqPostChangeNorm,[1 2]);
splits{6}       = strcmp(trialData.trialType,'Y') & trialData.audioFreqChangeNorm==2 & ismember(trialData.audioFreqPostChangeNorm,[3 4]);
splits{7}       = strcmp(trialData.trialType,'Y') & trialData.audioFreqChangeNorm==3 & ismember(trialData.audioFreqPostChangeNorm,[1 2]);
splits{8}       = strcmp(trialData.trialType,'Y') & trialData.audioFreqChangeNorm==3 & ismember(trialData.audioFreqPostChangeNorm,[3 4]);
nSplits         = length(splits);

% splits          = {};
% splits{1}       = strcmp(trialData.trialType,'X') & trialData.visualOriChangeNorm==3 & ismember(trialData.vecResponse,[1 2]);
% splits{2}       = strcmp(trialData.trialType,'X') & trialData.visualOriChangeNorm==3 & ~ismember(trialData.vecResponse,[1 2]);
% splits{3}       = strcmp(trialData.trialType,'Y') & trialData.audioFreqChangeNorm==3 & ismember(trialData.vecResponse,[1 2]);
% splits{4}       = strcmp(trialData.trialType,'Y') & trialData.audioFreqChangeNorm==3 & ~ismember(trialData.vecResponse,[1 2]);
% nSplits         = length(splits);
% splits{7}       = strcmp(trialData.trialType,'Y') & trialData.audioFreqChangeNorm==3 & ismember(trialData.audioFreqPostChangeNorm,[1 2]) & ~ismember(trialData.vecResponse,[1 2]);
% splits{8}       = strcmp(trialData.trialType,'Y') & trialData.audioFreqChangeNorm==3 & ismember(trialData.audioFreqPostChangeNorm,[3 4]) & ~ismember(trialData.vecResponse,[1 2]);

dims                = round(proc.ROI{1}{1}([3 4]));

output              = NaN(dims(2),dims(1),params.nTimebins_video,nSplits);
      
for iS = 1:nSplits
    for iT = 1:params.nTimebins_video
        %         selec       = squeeze(hist_mat(iT,:,splits{iS}));
        
%         selec       = squeeze(hits_mat_z_svd(iT,:,splits{iS}));
        
        selec       = squeeze(hist_mat(iT,:,splits{iS}));

        
        tempmot     = selec' * proc.uMotMask{1}'; %If A is an m-by-p and B is a p-by-n matrix, then C is an m-by-n matrix defined by
        %         tempmot     = selec' * (proc.uMotMask{1} - proc.avgmot{1})'; %If A is an m-by-p and B is a p-by-n matrix, then C is an m-by-n matrix defined by
        %         tempmot     = nanmean(abs(tempmot),1);
        %         tempmot     = nanmean(tempmot,1) - proc.avgmot{1}';
        tempmot        = reshape(tempmot,sum(splits{iS}),dims(2),dims(1));
        
%         if params.smoothing
%             for iTr = 1:sum(splits{iS})
%                 tempmot(iTr,:,:)          = conv2(squeeze(tempmot(iTr,:,:)),spat_filter,'same');
%             end
%         end
        temp     = squeeze(nanmean(tempmot,1));
        
        % dims = round(proc.ROI{1}{1}([3 4])) + 1;
        %         temp        = reshape(tempmot,dims(2),dims(1));
        
        output(:,:,iT,iS) = temp;
    end
end

%Make output absolute:
output = abs(output);

%Spatial filter:
spat_filter     = fspecial('gaussian',[10 10],params.spat_filter);

for iS = 1:nSplits
    for iT = 1:params.nTimebins_video
         output(:,:,iT,iS)          = conv2(output(:,:,iT,iS),spat_filter,'same');
    end
end

%upsample:
output2 = NaN(heightPixels,widthPixels,params.nTimebins_video,nSplits);
for iS = 1:nSplits
    for iT = 1:params.nTimebins_video
%          output2(:,:,iT,iS)          = interp2(output(:,:,iT,iS),linspace(1,dims(2),heightPixels),transpose(linspace(1,dims(1),widthPixels)));
         output2(:,:,iT,iS)          = interp2(output(:,:,iT,iS),transpose(linspace(1,dims(1),widthPixels)),linspace(1,dims(2),heightPixels));
    end
end

% figure; subplot(1,2,1); imagesc(output(:,:,iT,iS)); subplot(1,2,2); imagesc(output2(:,:,iT,iS));


%% Subtract baseline movement across video: 
% output2 = output2 - repmat(nanmean(nanmean(output2(:,:,params.xtime_video<0,:),4),3),1,1,params.nTimebins_video,nSplits);

bsl_mean    = repmat(nanmean(nanmean(output2(:,:,params.xtime_video<0,:),4),3),1,1,params.nTimebins_video,nSplits);
% temp        = output2(:,:,params.xtime_video<0,:);
% bsl_std     = repmat(nanstd(temp(:)),heightPixels,widthPixels,params.nTimebins_video,nSplits);
% output2     = output2 - bsl;
output2     = (output2 ./ bsl_mean)*100;

% output2 = output2 - repmat(nanmean(nanmean(output2(:,:,params.xtime_video<0,:),4),3),1,1,params.nTimebins_video,nSplits);
% output2 = output2 - repmat(nanmean(nanmean(output2(:,:,params.xtime_video<0,:),4),3),1,1,params.nTimebins_video,nSplits);

%% Figure:
iS = [5 6 7 8];
% iS = [7 8];

% clims = [prctile(output(:),params.clims(1)) prctile(output(:),params.clims(2))];

% clims = [prctile(output2(:),params.clims(1)) prctile(output2(:),params.clims(2))];
% clims = [130 prctile(output2(:),params.clims(2))];
% clims = [150 round(prctile(output2(:),99),-2)];
% clims = [150 600];
clims = params.clims;
idx = find(ismember(params.xtime_video,params.t_exframes));

cap.set(int8(1), params.framenumber-1);
%     cap.set(int8(1), iFr);
ret     = cap.read();
mt      = cell(ret(2));
x       = mt{1};
data    = double(py.array.array('d',py.numpy.nditer(x))); %d is for double, see link below on types
data    = data(1:3:end);
data    = reshape(data,[widthPixels heightPixels]);

cmap    = hot;

imagedat = [];

for iT = 1:length(idx)
    f = figure; hold all; set(gcf,'units','normalized','Position',[0.2 0.05 0.4 0.8],'color','w');
    figure(f)
%     B = nanmean(output(:,:,idx(iT),iS),4)';
    B = nanmean(output2(:,:,idx(iT),iS),4)';

    A = flipud(data);
    B = flipud(B);
    
    [hF,hB] = imoverlay(uint8(A),B,clims,[],'hot',0.6,f);
    colormap(cmap)
    
    title('')
    set(gca,'XTick',[],'YTick',[]) 
    xlabel(''); ylabel('');
    set(gca,'Visible','off')
    set(gca,'YDir','reverse')
    
    F = getframe(gca);
    imagedat(:,:,:,iT) = F.cdata;
    close(f)
end

f = figure; hold all; set(gcf,'units','normalized','Position',[0.2 0.4 0.1*length(idx) 0.39],'color','w');
for iT = 1:length(idx)
    ha = subplot(2,length(idx),iT); hold all;
    imagesc(uint8(imagedat(:,:,:,iT)));
%     imshow(imagedat(:,:,:,iT));
%     imshow(F.cdata);
    colormap(cmap)

    title(sprintf('%3.0f ms',params.xtime_video(idx(iT))*1e-3),'FontSize',10)
    set(gca,'XTick',[],'YTick',[]) 
    xlabel(''); ylabel('');
    set(gca,'Visible','off')
%     set(gca,'YDir','normal')
    set(gca,'YDir','reverse')
    
    c = colorbar();
    c.Ticks = [0 1];
    c.TickLabels = clims;
end

temp = squeeze(nanmean(nanmean(nanmean(output2(:,:,:,iS),4),2),1));

% if numel(iS)==1
%     temp = nanmean(hist_mat_z(:,splits{iS}),2);
% else
%     idx_temp = any([splits{iS}],2);
%     temp = nanmean(hist_mat_z(:,idx_temp),2);
% end

% f = figure; hold all; set(gcf,'units','normalized','Position',[0.1 0.4 0.32 0.18],'color','w');
subplot(2,length(idx),(1:length(idx))+length(idx)); hold all;
plot(params.xtime_video,temp,'k-','LineWidth',0.5)
plot(params.xtime_video(idx),temp(idx),'.','MarkerSize',25,'Color',[1 0 1])
plot([0 0],[0 ceil(max(temp))],':','Color',[0.6 0.6 0.6],'LineWidth',1)
xlim([params.t_pre params.t_post])
ylim([0 ceil(max(temp))])
sc = scalebar;
sc.XLen = 0.5e6;
sc.XUnit = 'sec';
sc.YLen = 100;
sc.YUnit = 'Change video ME (%)';
set(gca,'visible','off') %remove axes, only trace

filename = sprintf('MotionSVD_%s_%s.eps',params.session,params.experiment);
export_fig(fullfile(params.savedir,filename),gcf);

% export_fig(fullfile(params.savedir,filename),gcf,'-depsc');
% filename = sprintf('MotionSVD_%s_%s.eps',params.session,params.experiment);
% export_fig(fullfile(params.savedir,filename),gcf,-'opengl');


end


% %% Figure:
% f = figure; hold all; set(gcf,'units','normalized','Position',[0.2 0.4 0.32 0.3],'color','w');
% 
% iS = [5 6 7 8];
% 
% % clims = [prctile(output(:),params.clims(1)) prctile(output(:),params.clims(2))];
% 
% % clims = [prctile(output2(:),params.clims(1)) prctile(output2(:),params.clims(2))];
% % clims = [130 prctile(output2(:),params.clims(2))];
% clims = [150 round(prctile(output2(:),99),-2)];
% % clims = [150 600];
% % clims = params.clims;
% idx = find(ismember(params.xtime_video,params.t_exframes));
% 
% cap.set(int8(1), params.framenumber-1);
% %     cap.set(int8(1), iFr);
% ret     = cap.read();
% mt      = cell(ret(2));
% x       = mt{1};
% data    = double(py.array.array('d',py.numpy.nditer(x))); %d is for double, see link below on types
% data    = data(1:3:end);
% data    = reshape(data,[widthPixels heightPixels]);
% 
% cmap    = hot;
% 
% 
% for iT = 1:length(idx)
%     
%     ha = subplot(2,length(idx),iT); hold all;
%     
% %     B = nanmean(output(:,:,idx(iT),iS),4)';
%     B = nanmean(output2(:,:,idx(iT),iS),4)';
% 
%     A = flipud(data);
%     B = flipud(B);
%     
%     [hF,hB] = imoverlay(uint8(A),B,clims,[],'hot',0.6,ha);
%     colormap(cmap)
%     
%     title('')
%     set(gca,'XTick',[],'YTick',[]) 
%     xlabel(''); ylabel('');
%     set(gca,'Visible','off')
% %     set(gca,'YDir','normal')
%     set(gca,'YDir','reverse')
%           
%     F = getframe(ha);
%     cla;
%     imagesc(F.cdata);
% %     text(40,680,sprintf('%3.0f ms',params.xtime_video(idx(iT))*1e-3),'FontSize',12,'Color','w')
%     title(sprintf('%3.0f ms',params.xtime_video(idx(iT))*1e-3),'FontSize',10)
%     c = colorbar();
%     c.Ticks = clims;
% end
% 
% if numel(iS)==1
%     temp = nanmean(hist_mat_z(:,splits{iS}),2);
% else
%     idx_temp = any([splits{iS}],2);
%     temp = nanmean(hist_mat_z(:,idx_temp),2);
% end
% 
% % f = figure; hold all; set(gcf,'units','normalized','Position',[0.1 0.4 0.32 0.18],'color','w');
% ha = subplot(2,length(idx),(1:length(idx))+length(idx)); hold all;
% plot(params.xtime_video,temp,'k-','LineWidth',0.5)
% plot(params.xtime_video(idx),temp(idx),'.','MarkerSize',25,'Color',[1 0 1])
% plot([0 0],[0 0.5],':','Color',[0.6 0.6 0.6],'LineWidth',0.5)
% sc = scalebar;
% sc.XLen = 0.5e6;
% sc.XUnit = 'sec';
% sc.YLen = 0.5;
% sc.YUnit = 'z-score';
% set(gca,'visible','off') %remove axes, only trace
% 
% filename = sprintf('MotionSVD_%s_%s.eps',params.session,params.experiment);
% export_fig(fullfile(params.savedir,filename),gcf);
% % export_fig(fullfile(params.savedir,filename),gcf,'-depsc');
% % filename = sprintf('MotionSVD_%s_%s.eps',params.session,params.experiment);
% % export_fig(fullfile(params.savedir,filename),gcf,-'opengl');

