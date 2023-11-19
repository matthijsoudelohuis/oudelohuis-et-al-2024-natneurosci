function MOL_plotRasterVid_auV1_v2(sessionData,trialData,spikeData,video_zmat,cell_IDs,params) %#ok<INUSL>

nNeurons                = length(cell_IDs);

params.AlignOn          = 'stimChange';
params.tickresol        = 200e3;
params.sortBy           = 'motion';

% strings_cmaps = {'Blues' 'Blues' 'Reds' 'Reds'};

for iNeuron = 1:nNeurons %Loop over all neurons:
    
    cell_idx = strcmp(spikeData.cell_ID,cell_IDs(iNeuron));
    
    if any(cell_idx)
        %Get the relevant data for each session individually:
        temptrialData        = MOL_getTempPerSes(spikeData.session_ID(cell_idx),trialData);
        
        switch params.trialcategories
            case 'PostVolume'
                
                splits          = {};
                splits{1,1}       = ismember(temptrialData.trialType,{'X'}) & temptrialData.visualContrPostChangeNorm==1 & temptrialData.visualOriPostChangeNorm==1;
                splits{1,2}       = ismember(temptrialData.trialType,{'X'}) & temptrialData.visualContrPostChangeNorm==2 & temptrialData.visualOriPostChangeNorm==1;
                splits{1,3}       = ismember(temptrialData.trialType,{'X'}) & temptrialData.visualContrPostChangeNorm==2 & temptrialData.visualOriPostChangeNorm==2;
                splits{1,4}       = ismember(temptrialData.trialType,{'X'}) & temptrialData.visualContrPostChangeNorm==2 & temptrialData.visualOriPostChangeNorm==3;
                splits{1,5}       = ismember(temptrialData.trialType,{'X'}) & temptrialData.visualContrPostChangeNorm==3 & temptrialData.visualOriPostChangeNorm==1;
                
                splits{2,1}       = ismember(temptrialData.trialType,{'Y'}) & temptrialData.audioVolPostChangeNorm==1   & temptrialData.audioFreqPostChangeNorm==1;
                splits{2,2}       = ismember(temptrialData.trialType,{'Y'}) & temptrialData.audioVolPostChangeNorm==2   & temptrialData.audioFreqPostChangeNorm==1;
                splits{2,3}       = ismember(temptrialData.trialType,{'Y'}) & temptrialData.audioVolPostChangeNorm==2   & temptrialData.audioFreqPostChangeNorm==2;
                splits{2,4}       = ismember(temptrialData.trialType,{'Y'}) & temptrialData.audioVolPostChangeNorm==2   & temptrialData.audioFreqPostChangeNorm==3;
                splits{2,5}       = ismember(temptrialData.trialType,{'Y'}) & temptrialData.audioVolPostChangeNorm==3   & temptrialData.audioFreqPostChangeNorm==1;
                
                params.colors_splits       = {[163 12 235] [1 59 235] [0 25 107] [71 110 238] [0 235 195];     [235 143 12] [235 23 0] [193 0 235] [107 11 0] [238 88 71]};
                params.colors_splits       = cellfun(@(x) x/256,params.colors_splits,'UniformOutput',false);
                params.labels_splits       = {'VLow1' 'VMid1' 'VMid2' 'VMid3' 'VHigh1'; 'ALow1' 'AMid1' 'AMid2' 'AMid3' 'AHigh1' };
                params.lines_splits        = repmat({''},size(splits));
                
            case 'PostFeature'
                
                temptrialData.visualOriPostChangeNorm = round(temptrialData.visualOriPostChangeNorm/2);
                temptrialData.audioFreqPostChangeNorm = round(temptrialData.audioFreqPostChangeNorm/2);
                
                splits          = {};
                splits{1,1}       = ismember(temptrialData.trialType,{'X'}) & temptrialData.visualOriChangeNorm==2 & temptrialData.visualOriPostChangeNorm==1;
                splits{1,2}       = ismember(temptrialData.trialType,{'X'}) & temptrialData.visualOriChangeNorm==2 & temptrialData.visualOriPostChangeNorm==2;
                splits{1,3}       = ismember(temptrialData.trialType,{'X'}) & temptrialData.visualOriChangeNorm==3 & temptrialData.visualOriPostChangeNorm==1;
                splits{1,4}       = ismember(temptrialData.trialType,{'X'}) & temptrialData.visualOriChangeNorm==3 & temptrialData.visualOriPostChangeNorm==2;
                
                splits{2,1}       = ismember(temptrialData.trialType,{'Y'}) & temptrialData.audioFreqChangeNorm==2  & temptrialData.audioFreqPostChangeNorm==1;
                splits{2,2}       = ismember(temptrialData.trialType,{'Y'}) & temptrialData.audioFreqChangeNorm==2  & temptrialData.audioFreqPostChangeNorm==2;
                splits{2,3}       = ismember(temptrialData.trialType,{'Y'}) & temptrialData.audioFreqChangeNorm==3  & temptrialData.audioFreqPostChangeNorm==1;
                splits{2,4}       = ismember(temptrialData.trialType,{'Y'}) & temptrialData.audioFreqChangeNorm==3  & temptrialData.audioFreqPostChangeNorm==2;
                
                params.colors_splits       = {[12 220 200]  [30 110 255] [0 100 90] [20 0 107] ;     [235 143 12] [235 23 0] [193 0 235] [238 88 71]};
                params.colors_splits       = cellfun(@(x) x/256,params.colors_splits,'UniformOutput',false);
                
                params.labels_splits       = {'VthrAB' 'VthrCD' 'VmaxAB' 'VmaxCD'; 'AthrAB' 'AthrCD' 'AmaxAB' 'AmaxCD';};
                params.lines_splits        = repmat({''},size(splits));
                
            case 'PostFeature2'
                
                temptrialData.visualOriPostChangeNorm = round(temptrialData.visualOriPostChangeNorm/2);
                temptrialData.audioFreqPostChangeNorm = round(temptrialData.audioFreqPostChangeNorm/2);
                
                splits          = {};
                splits{1,1}       = ismember(temptrialData.trialType,{'X'}) & temptrialData.visualOriChangeNorm==3  & temptrialData.visualOriPostChangeNorm==1;
                splits{1,2}       = ismember(temptrialData.trialType,{'X'}) & temptrialData.visualOriChangeNorm==3  & temptrialData.visualOriPostChangeNorm==2;
                
                splits{2,1}       = ismember(temptrialData.trialType,{'Y'}) & temptrialData.audioFreqChangeNorm==3   & temptrialData.audioFreqPostChangeNorm==1;
                splits{2,2}       = ismember(temptrialData.trialType,{'Y'}) & temptrialData.audioFreqChangeNorm==3   & temptrialData.audioFreqPostChangeNorm==2;
                
                params.colors_splits       = {[0 100 90] [20 0 107] ;     [193 0 235] [238 88 71]};
                params.colors_splits       = cellfun(@(x) x/256,params.colors_splits,'UniformOutput',false);
                
%                 params.labels_splits       = {'VthrAB' 'VthrCD' 'VmaxAB' 'VmaxCD'; 'AthrAB' 'AthrCD' 'AmaxAB' 'AmaxCD';};
                params.labels_splits       = {'VmaxAB' 'VmaxCD'; 'AmaxAB' 'AmaxCD';};
                params.lines_splits        = repmat({''},size(splits));
                
            case 'Saliency'
                                
                splits          = {};
                splits{1,1}       = ismember(temptrialData.trialType,{'X'}) & temptrialData.visualOriChangeNorm==2;
                splits{1,2}       = ismember(temptrialData.trialType,{'X'}) & temptrialData.visualOriChangeNorm==3;
                
                splits{2,1}       = ismember(temptrialData.trialType,{'Y'}) & temptrialData.audioFreqChangeNorm==2;
                splits{2,2}       = ismember(temptrialData.trialType,{'Y'}) & temptrialData.audioFreqChangeNorm==3;
                
%                 params.colors_splits       = {[12 220 200]  [30 110 255] [0 100 90] [20 0 107];     [235 143 12] [235 23 0] [193 0 235] [238 88 71]};
                params.colors_splits       = {[12 220 200]  [30 110 255];     [235 143 12] [235 23 0]};
                params.colors_splits       = cellfun(@(x) x/256,params.colors_splits,'UniformOutput',false);
                
                params.labels_splits       = {'Vthr' 'Vmax'; 'Athr' 'Amax';};
                params.lines_splits        = repmat({''},size(splits));
                
            case 'State'
                
%                 splits          = {};
%                 splits{1,1}       = ismember(temptrialData.trialType,{'Y'}) & temptrialData.visualOriChangeNorm==2;
%                 splits{1,2}       = ismember(temptrialData.trialType,{'Y'}) & temptrialData.visualOriChangeNorm==3;
%                 
%                 splits{2,1}       = ismember(temptrialData.trialType,{'Y'}) & temptrialData.audioFreqChangeNorm==2;
%                 splits{2,2}       = ismember(temptrialData.trialType,{'Y'}) & temptrialData.audioFreqChangeNorm==3;
                
                splits          = {};
%                 splits{1,1}       = ismember(temptrialData.trialType,{'X'}) & temptrialData.visualOriChangeNorm==2 & temptrialData.State==1;
%                 splits{1,2}       = ismember(temptrialData.trialType,{'X'}) & temptrialData.visualOriChangeNorm==3 & temptrialData.State==1;
                
%                 splits{1,1}       = ismember(temptrialData.trialType,{'Y'}) & temptrialData.audioFreqChangeNorm==2 & temptrialData.State==1;
%                 splits{1,2}       = ismember(temptrialData.trialType,{'Y'}) & temptrialData.audioFreqChangeNorm==3 & temptrialData.State==1;
%                 
%                 splits{2,1}       = ismember(temptrialData.trialType,{'Y'}) & temptrialData.audioFreqChangeNorm==2 & temptrialData.State==2;
%                 splits{2,2}       = ismember(temptrialData.trialType,{'Y'}) & temptrialData.audioFreqChangeNorm==3 & temptrialData.State==2;
                
                
                splits{1,1}       = ismember(temptrialData.trialType,{'X'}) & temptrialData.State==1;
                splits{1,2}       = ismember(temptrialData.trialType,{'Y'}) & temptrialData.State==1;
                
                splits{2,1}       = ismember(temptrialData.trialType,{'X'}) & temptrialData.State==2;
                splits{2,2}       = ismember(temptrialData.trialType,{'Y'}) & temptrialData.State==2;
                
%                 params.colors_splits       = {[12 220 200]  [30 110 255] [0 100 90] [20 0 107];     [235 143 12] [235 23 0] [193 0 235] [238 88 71]};
%                 params.colors_splits       = {[12 220 200]  [30 110 255];     [235 143 12] [235 23 0]};
                params.colors_splits       = {[235 143 12] [235 23 0];     [235 143 12] [235 23 0]};
                params.colors_splits       = cellfun(@(x) x/256,params.colors_splits,'UniformOutput',false);
                
                params.labels_splits       = {'Athr' 'Amax'; 'Athr' 'Amax';};
                params.lines_splits        = repmat({''},size(splits));
                
        end
        
        if all(cellfun(@sum,splits)>=params.minNtrials,[1 2])
            [params.nSplitsX,params.nSplitsY] = size(splits);
            
            %Get times:
            events_ts               = temptrialData.(params.AlignOn);
            spikes_ts               = spikeData.ts{cell_idx};
            
            %sort spikes
            vecTrialPerSpike        = nan(size(spikes_ts));
            vecTimePerSpike         = nan(size(spikes_ts));
            
            % Get spiking times:
            for ev = 1:length(events_ts) % Get spikes within window around event ev:
                idx                     = spikes_ts>events_ts(ev)+params.t_pre & spikes_ts<=events_ts(ev)+params.t_post;
                vecTrialPerSpike(idx)   = ev;
                vecTimePerSpike(idx)    = spikes_ts(idx)-events_ts(ev);
            end
            
            %Compute histogram:
            events_ts               = temptrialData.(params.AlignOn);
            hist_mat                = calc_psth(events_ts,spikes_ts,params);    %Construct histogram matrix
            
            ratemat                 = NaN(params.nSplitsX,params.nSplitsY,params.nTimebins); %compute average per condition
            for iSx = 1:params.nSplitsX %Store the mean response for each of these splits
                for iSy = 1:params.nSplitsY %Store the mean response for each of these splits
                    ratemat(iSx,iSy,:) = nanmean(hist_mat(splits{iSx,iSy},:),1); %#ok<*AGROW>
                end
            end
            
            ratemax = ceil(max(ratemat(:)));
            
            %Video:
            idx_ses                 = strcmp(trialData.session_ID,spikeData.session_ID(cell_idx));
            movz                    = video_zmat(:,idx_ses)';
            
            nTrials                 = size(movz,1);
            if params.smoothSVD
                for iTrial = 1:nTrials %some smoothing for visualization
                    %                 temp_movz(iTrial,:) = smooth(temp_movz(iTrial,:),params.smoothSVD,'sgolay');
                    movz(iTrial,:) = smooth(movz(iTrial,:),params.smoothSVD,'sgolay');
                end
            end
            
            motmat                 = NaN(params.nSplitsX,params.nSplitsY,params.nTimebins_video); %compute average per condition
            for iSx = 1:params.nSplitsX %Store the mean response for each of these splits
                for iSy = 1:params.nSplitsY %Store the mean response for each of these splits
                    motmat(iSx,iSy,:) = nanmean(movz(splits{iSx,iSy},:),1); %#ok<*AGROW>
                end
            end
            
            motlim = [-1 prctile(movz(:),97)];
            
            figure; set(gcf,'units','normalized','Position',[0.1 0.4 0.38 0.37],'color','w'); hold all;
            xwidth = 0.2;
            subplot(2,params.nSplitsX*2,1); hold all;
            set(gca,'Position',[0.03 0.7 xwidth 0.25])
            handles = [];
            iSx = 1; %Store the mean response for each of these splits
            for iSy = 1:params.nSplitsY %Store the mean response for each of these splits
                handles(iSy) = plot(params.xtime,squeeze(ratemat(iSx,iSy,:)),'-','Color',params.colors_splits{iSx,iSy},'LineWidth',0.5);
            end
            
            xlabel('')
            xlim([params.t_pre+100e3 params.t_post]);
            ylim([0 ratemax])
            set(gca,'XTick',params.t_pre:params.tickresol:params.t_post,'XTickLabels',[],'YTick',get(gca,'YLim'),'FontSize',8)
            set(gca,'XMinorTick','on','TickDir', 'out')
            plot([0 0],[0 ratemax],'k--','LineWidth',0.5)
            ylabel('')
            title(cell_IDs(iNeuron),'FontSize',10)
            legend(handles,params.labels_splits(iSx,:),'Location','NorthEast','FontSize',8); legend boxoff;
            
            subplot(2,params.nSplitsX*2,2); hold all;
            set(gca,'Position',[0.03+xwidth*1.1 0.7 xwidth 0.25])
            handles = [];
            iSx = 1; %Store the mean response for each of these splits
            for iSy = 1:params.nSplitsY %Store the mean response for each of these splits
                handles(iSy) = plot(params.xtime_video,squeeze(motmat(iSx,iSy,:)),'-','Color',params.colors_splits{iSx,iSy},'LineWidth',0.5);
            end
            
            xlabel('')
            xlim([params.t_pre+100e3 params.t_post]);
            ylim(motlim)
            set(gca,'XTick',params.t_pre:params.tickresol:params.t_post,'XTickLabels',[],'YTick',get(gca,'YLim'),'FontSize',8)
            set(gca,'XMinorTick','on','TickDir', 'out')
            plot([0 0],[-1 ratemax],'k--','LineWidth',0.5)
            ylabel('')
            legend(handles,params.labels_splits(iSx,:),'Location','NorthEast','FontSize',8); legend boxoff;
            
            subplot(2,params.nSplitsX*2,3); hold all;
            set(gca,'Position',[0.03+xwidth*1.1*2 0.7 xwidth 0.25])
            handles = [];
            iSx = 2; %Store the mean response for each of these splits
            for iSy = 1:params.nSplitsY %Store the mean response for each of these splits
                handles(iSy) = plot(params.xtime,squeeze(ratemat(iSx,iSy,:)),'-','Color',params.colors_splits{iSx,iSy},'LineWidth',0.5);
            end
            
            xlabel('')
            xlim([params.t_pre+100e3 params.t_post]);
            ylim([0 ratemax])
            set(gca,'XTick',params.t_pre:params.tickresol:params.t_post,'XTickLabels',[],'YTick',get(gca,'YLim'),'FontSize',8)
            set(gca,'XMinorTick','on','TickDir', 'out')
            plot([0 0],[0 ratemax],'k--','LineWidth',0.5)
            ylabel('')
            legend(handles,params.labels_splits(iSx,:),'Location','NorthEast','FontSize',8); legend boxoff;
            
            subplot(2,params.nSplitsX*2,4); hold all;
            set(gca,'Position',[0.03+xwidth*1.1*3 0.7 xwidth 0.25])
            handles = [];
            iSx = 2; %Store the mean response for each of these splits
            for iSy = 1:params.nSplitsY %Store the mean response for each of these splits
                handles(iSy) = plot(params.xtime_video,squeeze(motmat(iSx,iSy,:)),'-','Color',params.colors_splits{iSx,iSy},'LineWidth',0.5);
            end
            
            xlabel('')
            xlim([params.t_pre+100e3 params.t_post]);
            ylim([-1 4])
            set(gca,'XTick',params.t_pre:params.tickresol:params.t_post,'XTickLabels',[],'YTick',get(gca,'YLim'),'FontSize',8)
            set(gca,'XMinorTick','on','TickDir', 'out')
            plot([0 0],[-1 ratemax],'k--','LineWidth',0.5)
            ylabel('')
            legend(handles,params.labels_splits(iSx,:),'Location','NorthEast','FontSize',8); legend boxoff;
            
            for iSx = 1:params.nSplitsX
                
                idx_trials = [];
                idx_trials_cell = {};
                for iSy = 1:params.nSplitsY %Store the mean response for each of these splits
                    
                    %Sort trials based on post-change orientation and response latency (orientation is sorted by adding a lot of latency)
                    trialselec              = find(splits{iSx,iSy}');
                    switch params.sortBy
                        case 'responseLatency'
                            temp                    = temptrialData.responseLatency(trialselec);
                            temp(isnan(temp))       = 0;
                            [~,idx]                 = sort(temp,'descend');
                        case 'motion'
                            temp                    = nanmean(movz(trialselec,params.xtime_video>0 & params.xtime_video<1e6),2);
                            [~,idx]                 = sort(temp,'descend');
                    end
                    
                    idx_trials              = [idx_trials; trialselec(idx)'];
                    idx_trials_cell{iSy}    = trialselec(idx)';
                end
                
                %Plot heatmap of firing rate:
                hSub = subplot(2,params.nSplitsX*2,4+iSx*2-1); hold all;
                set(gca,'Position',[0.03+xwidth*1.1*(iSx-1)*2 0.1 xwidth 0.56])
                
                temp_rate   = flipud(hist_mat(idx_trials,:))';
                
                [X,Y]       = meshgrid(params.xtime,1:length(idx_trials));
                h           = pcolor(X,Y,temp_rate');
                set(h, 'EdgeColor', 'none');
                
                if iSx==2
                    colormap(hSub,getPyPlot_cMap('Reds'));
                else
                    colormap(hSub,getPyPlot_cMap('Blues'));
                end
                
                trialcounter = length(idx_trials);
                for iSy = 1:params.nSplitsY 
                    nTrials             = length(idx_trials_cell{iSy});
                    plot([params.t_pre+100e3 params.t_pre+100e3],[trialcounter+1 trialcounter-nTrials+1],'-','Color',params.colors_splits{iSx,iSy},'LineWidth',2)
                    trialcounter        = trialcounter-nTrials;
                end
                
                %                     temptc = trialcounter;
                %                     for intTrial = 1:length(idx_trials_cell{iSy})
                %                         lat         = temptrialData.responseLatency(idx_trials_cell{iSy}(intTrial));
                %                         if ~isnan(lat)
                % %                             line([lat lat],[temptc;temptc-1],'Color','k','LineWidth',2);
                %                         end
                %                         temptc = temptc-1;
                %                     end
                %                     trialcounter        = trialcounter-nTrials;
                %
                
                lims        = [0 prctile(hist_mat(:),97)];
                caxis(lims)
                plot([0 0],[0 length(idx_trials)],'k--','LineWidth',0.5)
                ylim([1 length(idx_trials)]);
                xlim([params.t_pre+100e3 params.t_post]);
                xlabel('Time (ms)','FontSize',9);
                set(gca,'XTick',params.t_pre:params.tickresol:params.t_post,'XTickLabels',(params.t_pre:params.tickresol:params.t_post)*1e-3,'YTick',[],'FontSize',8)
                set(gca,'XMinorTick','on','TickDir', 'out')
                
                %Video motion:
                hSub = subplot(2,params.nSplitsX*2,4+iSx*2); hold all;
                set(gca,'Position',[0.03+xwidth*1.1*((iSx-1)*2+1) 0.1 xwidth 0.56])
                
                temp_movz   = flipud(movz(idx_trials,:))';
                
                [X,Y]       = meshgrid(params.xtime_video,1:length(idx_trials));
                h           = pcolor(X,Y,temp_movz');
                set(h, 'EdgeColor', 'none');
                
%                 trialcounter = length(idx_trials);
%                 for iSy = 1:params.nSplitsY 
%                     trialselec          = find(splits{iSx,iSy}');
%                     temp                = temptrialData.responseLatency(trialselec);
%                     temp(isnan(temp))   = 0;
%                     orig                = trialcounter;%+length(trialselec);
% %                     patch([params.t_post params.t_post params.t_post-12.5e3 params.t_post-12.5e3],[orig orig-sum(temp~=0) orig-sum(temp~=0) orig],[0 0 0],...
% %                         'EdgeColor',[0.6 0.6 0.6],'FaceColor',[0.6 0.6 0.6],'FaceAlpha',1,'EdgeAlpha',1)
% 
%                     trialcounter        = trialcounter-length(trialselec);
%                 end
                
                colormap(hSub,getPyPlot_cMap('CMRmap'));
                
                caxis(motlim)
                plot([0 0],[0 length(idx_trials)],'k--','LineWidth',0.5)
                ylim([1 length(idx_trials)]);
                xlim([params.t_pre+100e3 params.t_post]);
                xlabel('Time (ms)','FontSize',9);
                set(gca,'XTick',params.t_pre:params.tickresol:params.t_post,'XTickLabels',(params.t_pre:params.tickresol:params.t_post)*1e-3,'YTick',[],'FontSize',8)
                set(gca,'XMinorTick','on','TickDir', 'out')
            end
            
            if isfield(params,'exportfig') && params.exportfig
                export_fig(fullfile(params.savedir,sprintf('ExNeuron_%s',spikeData.cell_ID{cell_idx})),'-eps','-nocrop')
                %                 export_fig(fullfile(params.savedir,sprintf('ExNeuron_%s',spikeData.cell_ID{cell_idx})),'-png','-nocrop')
            end
        end
        
    end
    
end


end