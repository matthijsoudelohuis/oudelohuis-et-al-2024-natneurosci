
function plotDominanceMap(DominanceIndexMat,DIlims)

% DominanceIndexMat   = flipud(DominanceIndexMat);

nCond               = size(DominanceIndexMat,1);
params.colormap     = 'redblue';

figure; set(gcf,'color','w','units','normalized','Position', [0.1 0.5 .17 .2]);

switch nCond
    case 2
        axislabels = {'Catch' 'Thr'};
    case 3
        axislabels = {'Catch' 'Thr' 'Max'};
    case 5
        axislabels = {'Catch' 'Sub' 'Thr' 'Sup' 'Max'};
end

%% Show Heatmap with dominance indices:
% imagesc(DominanceIndexMat,[-1 1]); hold on;
pcolor([DominanceIndexMat nan(nCond,1); nan(1,nCond+1)]); %Show Heatmap with dominance indices
params.colormap     = 'redblue';

switch params.colormap
    case 'redblue'
        h       = coolwarm(64) + repmat(0.1-abs(linspace(-0.1,.1,64))',1,3);
        h(h>1)  = 1;
        colormap(h);
    case 'parula'
        colormap(parula);
end
set(gca,'XTick',1.5:nCond+0.5,'XTickLabels',axislabels,'XTickLabelRotation',45,'YTick',1.5:nCond+0.5,'YTickLabels',axislabels,'FontSize',10)
xlabel('Visual','FontSize',10)
ylabel('Auditory','FontSize',10)

c = colorbar('YTick',[DIlims(1) 0 DIlims(2)]);
c.Label.String = 'Auditory Dominance Index';
caxis(DIlims)

end