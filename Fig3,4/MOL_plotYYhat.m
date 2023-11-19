function MOL_plotYYhat(temptrialData,params,Yvars,cell_ID)

nYvars              = length(Yvars);

nTrials             = length(temptrialData.session_ID);

clims               = [2 prctile(Yvars{1}(:),97)];
% clims               = [2 prctile(Yvars{1}(:),99)];
% cmapdata            = parula(128);
cmapdata            = jet(128); cmapdata = cmapdata(1:end-5,:).^0.5;
% cmapdata            = turbo(128); cmapdata = cmapdata(1:end-5,:).^2;


% cmapdata            = getPyPlot_cMap('gist_rainbow');
% cmapdata            = getPyPlot_cMap('hot');
% cmapdata            = getPyPlot_cMap('jet');
% cmapdata            = getPyPlot_cMap('YlGnBu');
% cmapdata            = getPyPlot_cMap('rainbow');
% cmapdata            = getPyPlot_cMap('gist_stern');
% cmapdata = flipud(cmapdata);

temptrialData.visualOriPostChangeNorm = round(temptrialData.visualOriPostChangeNorm/2);
temptrialData.audioFreqPostChangeNorm = round(temptrialData.audioFreqPostChangeNorm/2);

splits          = {};
splits{1}       = ismember(temptrialData.trialType,{'X'}) & temptrialData.visualOriChangeNorm==2 & temptrialData.visualOriPostChangeNorm==1;
splits{2}       = ismember(temptrialData.trialType,{'X'}) & temptrialData.visualOriChangeNorm==2 & temptrialData.visualOriPostChangeNorm==2;
splits{3}       = ismember(temptrialData.trialType,{'X'}) & temptrialData.visualOriChangeNorm==3 & temptrialData.visualOriPostChangeNorm==1;
splits{4}       = ismember(temptrialData.trialType,{'X'}) & temptrialData.visualOriChangeNorm==3 & temptrialData.visualOriPostChangeNorm==2;

splits{5}       = ismember(temptrialData.trialType,{'Y'}) & temptrialData.audioFreqChangeNorm==2 & temptrialData.audioFreqPostChangeNorm==1;
splits{6}       = ismember(temptrialData.trialType,{'Y'}) & temptrialData.audioFreqChangeNorm==2 & temptrialData.audioFreqPostChangeNorm==2;
splits{7}       = ismember(temptrialData.trialType,{'Y'}) & temptrialData.audioFreqChangeNorm==3 & temptrialData.audioFreqPostChangeNorm==1;
splits{8}       = ismember(temptrialData.trialType,{'Y'}) & temptrialData.audioFreqChangeNorm==3 & temptrialData.audioFreqPostChangeNorm==2;

splits{9}       = ismember(temptrialData.trialType,{'P'});
% splits{9}       = ismember(temptrialData.trialType,{'P'}) & ismember(temptrialData.vecResponse,3);
% splits{10}      = ismember(temptrialData.trialType,{'P'}) & ismember(temptrialData.vecResponse,[1 2]);

nSplits         = length(splits);

% params.labels_splits       = {'VthrAB' 'VthrCD' 'VmaxAB' 'VmaxCD' 'AthrAB' 'AthrCD' 'AmaxAB' 'AmaxCD';};
params.labels_splits       = {'VthrAB' 'VthrCD' 'VmaxAB' 'VmaxCD' 'AthrAB' 'AthrCD' 'AmaxAB' 'AmaxCD' 'Catch';};
params.lines_splits        = repmat({''},size(splits));

figure; set(gcf,'units','normalized','Position',[0.2 0.13 0.6 0.73],'color','w'); hold all;

for iY = 1:nYvars
    for iS = 1:nSplits
        nSelecTrials    = sum(splits{iS});
        Y_r             = reshape(Yvars{iY},params.nTimebins,nTrials);
        Y_r             = Y_r(:,splits{iS});
        %total mode:
        subplot(nSplits,nYvars,iY + (iS-1)*nYvars); hold all;
        [Xax,Yax] = meshgrid(params.xtime,1:nSelecTrials);
        h           = pcolor(Xax,Yax,Y_r');
        set(h, 'EdgeColor', 'none');
        set(gca,'XTick',[],'YTick',[]);
        xlim([params.t_pre params.t_post])
        ylim([0 nSelecTrials])
        plot([0 0],[0 nSelecTrials],'w--','LineWidth',0.5)
        caxis(clims);
        colormap(cmapdata);
        if iS == 1
            title(params.Yvarlabels{iY})
        end
        if iY == 1
            ylabel(params.labels_splits{iS},'FontSize',9)
        end
        
        if iS == nSplits && iY == nYvars
            xlabel(cell_ID,'FontSize',9)
        end
        
    end
end

postemp = get(gca,'Position');
hb = colorbar('location','eastoutside');
set(gca,'position',postemp);


if isfield(params,'exportfig') && params.exportfig
    export_fig(fullfile(params.savedir,sprintf('ExNeuron_GLM_%s',cell_ID)),'-eps','-nocrop')
    %                 export_fig(fullfile(params.savedir,sprintf('ExNeuron_%s',spikeData.cell_ID{cell_idx})),'-png','-nocrop')
end


end