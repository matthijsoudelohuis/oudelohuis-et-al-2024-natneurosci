%% Oude Lohuis et al. 2024 Nat Neurosci
% Triple dissociation of auditory, visual, and motor processing in primary visual cortex
% MOL (C) 2023
% Show video motion across the face as energy

startover
load('Dataset2_2.mat')

%% 
figure; set(gcf,'units','normalized','Position',[0.45 0.4 0.24 0.46],'color','w')

temp = reshape(proc.avgframe,size(proc.wpix{1},1),size(proc.wpix{1},2));

imagesc(flipud(temp'))

h       =  coolwarm(256)  + repmat(0.1-abs(linspace(-0.1,.1,256))',1,3);
h(h>1)  = 1;
colormap(h);
ylabel(''); xlabel(''); set(gca,'XTick',[],'YTick',[]);

%% 
figure; set(gcf,'units','normalized','Position',[0.45 0.4 0.24 0.46],'color','w')

temp = reshape(proc.avgmotion,size(proc.wpix{1},1),size(proc.wpix{1},2));

imagesc(flipud(temp'))

h       =  coolwarm(256)  + repmat(0.1-abs(linspace(-0.1,.1,256))',1,3);
h(h>1)  = 1;
colormap(h);
ylabel(''); xlabel(''); set(gca,'XTick',[],'YTick',[]);

%% 
dims = round(proc.ROI{1}{1}([3 4])) + 1;
dims = round(proc.ROI{1}{1}([3 4]));
nSVDs = 3;

figure; set(gcf,'units','normalized','Position',[0.05 0.4 0.58 0.28],'color','w')
for iSVD = 1:nSVDs
    subplot(1,nSVDs,iSVD)
    temp = reshape(proc.uMotMask{1}(:,iSVD),dims(2),dims(1));
    imagesc(flipud(temp'))
    
    h       = coolwarm(256) + repmat(0.1-abs(linspace(-0.1,.1,256))',1,3);
    h(h>1)  = 1;
    colormap(h);
    ylabel(''); xlabel(''); set(gca,'XTick',[],'YTick',[]);
    title(sprintf('PCA %d',iSVD));
end


%%
figure; set(gcf,'units','normalized','Position',[0.3 0.15 0.6 0.25],'color','w'); hold all;

selec = 1:1000;

for iSVD = 1:nSVDs
    plot(proc.motSVD{1}(selec,iSVD),'LineWidth',2);
end


