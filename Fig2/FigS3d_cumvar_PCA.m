%% Script that analyzes sensory responses to auditory stimuli across passive and active task versions:
% as reported in Oude Lohuis et al. 2023 Nat Neurosci
% MOL (C) 2023

%%
startover
videoFiles = {};

%% Vids from separate animals:
videoFiles{end+1}       = 'F:\Data\CHDET\RawData\2028\2019-12-13_12-46-26\VisOnlyTwolevels\VT1.0001';
videoFiles{end+1}       = 'F:\Data\CHDET\RawData\2035\2019-12-16_15-38-43\VisOnlyTwolevels\VT1.0001';
videoFiles{end+1}       = 'F:\Data\CHDET\RawData\2044\2021-04-28_14-00-10\ChangeDetectionConflict\VT1.0001';
videoFiles{end+1}       = 'F:\Data\CHDET\RawData\2030\2020-01-21_14-41-38\ChangeDetectionConflict\VT1.0002';
videoFiles{end+1}       = 'F:\Data\CHDET\RawData\2029\2019-12-12_13-45-39\VisOnlyTwolevels\VT1.0002';
videoFiles{end+1}       = 'F:\Data\CHDET\RawData\2034\2019-12-18_13-31-12\VisOnlyTwolevels\VT1.0002';
videoFiles{end+1}       = 'F:\Data\CHDET\RawData\2045\2021-04-29_15-55-10\ChangeDetectionConflict\VT1.0001';
videoFiles{end+1}       = 'F:\Data\CHDET\RawData\2044\2021-04-23_13-12-02\ChangeDetectionConflict\VT1.0001';

%% Loop over data
nVids       = length(videoFiles);
g           = NaN(nVids,500);
for iVid = 1:nVids
    g(iVid,:) = compVarExpl(videoFiles{iVid});
end

%% Save 
save('Dataset2_4.mat',g)

%% or load preprocessed data:
load('Dataset2_4.mat')

%% Make figure and determine n video PCs based on threshold on var explained:
thresh = 0.8;
yticks = [0 .2 0.4 0.6 0.8 1];
g_norm = g ./ repmat(max(g,[],2),1,500);
fprintf('#PCs for %2.2f of normalized variance: \n%d\n',thresh,find(mean(g_norm,1)>thresh,1))
nCompUse = find(mean(g_norm,1)>thresh,1);
figure; set(gcf,'units','normalized','Position',[0.45 0.4 0.24 0.18],'color','w'); hold all;
subplot(1,2,1); hold all;

handles = [];
h = shadedErrorBar(1:500,mean(g,1),std(g,[],1),{'Color','r','LineWidth',1},0);
handles(1) = h.mainLine; delete(h.edge(1)); delete(h.edge(2));
h = shadedErrorBar(1:500,mean(g_norm,1),std(g_norm,[],1),{'Color','b','LineWidth',1},0);
handles(2) = h.mainLine; delete(h.edge(1)); delete(h.edge(2));

% plot(1:500,mean(g,1),'r-');
% plot(1:500,mean(g_norm,1),'b-');
plot([0 500],[thresh thresh],'k:')
plot([nCompUse nCompUse],[0 100],'b:')

xlim([0 500])
ylim([0 1])
% plot([25 25],[0 100],'k:')
set(gca,'XTick',[1 50 100 250 500] ,'YTick',yticks);
ylabel('Explained variance')
xlabel('#Principle Components')
legend(handles,{'Raw' 'Normalized'}); legend('Location','South'); legend boxoff;
subplot(1,2,2); hold all;

h = shadedErrorBar(1:500,mean(g,1),std(g,[],1),{'Color','r','LineWidth',1},0);
handles(1) = h.mainLine; delete(h.edge(1)); delete(h.edge(2));
h = shadedErrorBar(1:500,mean(g_norm,1),std(g_norm,[],1),{'Color','b','LineWidth',1},0);
handles(2) = h.mainLine; delete(h.edge(1)); delete(h.edge(2));

% plot(1:500,mean(g,1),'r-');
% plot(1:500,mean(g_norm,1),'b-');

xlim([0.5 50])
ylim([0 1])
plot([0.001 500],[thresh thresh],'k:')
plot([nCompUse nCompUse],[0 100],'b:')
set(gca,'XTick',sort([1 2 5 10 20 50 100 250 500 nCompUse]),'YTick',yticks);
set(gca,'Xscale','log')
fprintf('EV at 500 video PCs: \n%2.2f %%\n',mean(g(:,500),1)*100)
