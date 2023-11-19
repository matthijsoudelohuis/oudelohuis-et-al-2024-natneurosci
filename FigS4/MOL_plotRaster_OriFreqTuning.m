function MOL_plotRaster_OriFreqTuning(sessionData,trialData,spikeData,cell_IDs,params)

nNeurons                = length(cell_IDs);

params.nSplits          = 8;

for iNeuron = 1:nNeurons %Loop over all neurons:
    
    cell_idx = strcmp(spikeData.cell_ID,cell_IDs(iNeuron));
    
    if any(cell_idx)
        %Get the relevant data for each session individually:
        temptrialData        = MOL_getTempPerSes(spikeData.session_ID(cell_idx),trialData);
        
        params.trialcategories = 'OriFreq';
        
        switch params.trialcategories
            case 'OriFreq'
                splits          = {};
                for iOri = 1:params.nSplits
                    splits{iOri,1} = strcmp(temptrialData.trialType,'V') & temptrialData.visualOri==params.allOris(iOri);
                end
                
                for iFreq = 1:params.nSplits
                    splits{iFreq,2} = strcmp(temptrialData.trialType,'A') & temptrialData.audioFreq==params.allFreqs(iFreq);
                end
        end
        
        if all(all(cellfun(@sum,splits)>=params.minNtrials))
            
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
            hist_mat                = calc_psth(events_ts,spikes_ts,params);    %Construct histogram matrix:
            
            figure; set(gcf,'units','normalized','Position',[0.1 0.3 0.2 0.5],'color','w'); hold all;
            
            for iMod = 1:2
                ratemat = [];
                for iSplit = 1:params.nSplits %Store the mean response for each of these splits
                    ratemat(iSplit,:) = nanmean(hist_mat(splits{iSplit,iMod},:),1); %#ok<*AGROW>
                end
                
                subplot(params.nSplits*2+2,2,iMod); hold all;
                for iSplit = 1:params.nSplits %Store the mean response for each of these splits
%                     handles(iSplit) = plot(params.xtime,ratemat(iSplit,:),params.lines_ztrials{iSplit},'Color',params.colors_ztrials{iSplit},'LineWidth',3);
                    handles(iSplit) = plot(params.xtime,ratemat(iSplit,:),'- ','Color',params.colors_ztrials(iSplit + (iMod-1)*params.nSplits,:),'LineWidth',3);
                end
                if iMod==1
                    title(cell_IDs(iNeuron),'FontSize',12)
                end
                
                xlabel('')
                xlim([params.t_pre params.t_post]);
                ylim([0 ceil(max(max(ratemat)))])
                set(gca,'XTick',[],'YTick',get(gca,'YLim'))
                plot([0 0],[0 max(max(ratemat))],'k--','LineWidth',2)
%                 if iNeuron ==1
%                     ylabel('Firing rate (sp s-1)','FontSize', 15)
%                     legend(handles,params.labels_ztrials,'FontSize',10); legend boxoff
%                 else ylabel('')
%                 end
                 ylabel('')

                 for iSplit = 1:params.nSplits
                    subplot(params.nSplits+2,2,(iSplit*2)+iMod); hold all;
                    
                    %Sort trials based on post-change orientation and response latency (orientation is sorted by adding a lot of latency)
                    trialselec          = find(splits{iSplit,iMod}');
                    temp                = temptrialData.responseLatency(trialselec);
                    temp(isnan(temp))   = 0;
                    [~,idx]             = sort(temp,'descend');
                    trialselec          = trialselec(idx);
                    
                    %plot per trial
                    %             h = patch([params.t_pre params.t_post params.t_post params.t_pre],[0 0 length(trialselec) length(trialselec)],'k');
                    %             h.FaceAlpha = 0.1;
                    %             h.FaceAlpha = 1;
                    %             h.EdgeColor = [1 1 1];
                    
                    for intTrial = 1:length(trialselec)%find(trialselec') %1:numel(vecTrialStarts)
                        
                        vecTimes = vecTimePerSpike(vecTrialPerSpike==trialselec(intTrial));
                        lat = temptrialData.responseLatency(trialselec(intTrial));
                        if ~isnan(lat)
                            line([lat lat],[intTrial-0.5;intTrial+0.5],'Color','k','LineWidth',1);
                        end
                        line([vecTimes(:)';vecTimes(:)'],[intTrial*ones(1,numel(vecTimes))-1;intTrial*ones(1,numel(vecTimes))+1],'Color',params.colors_ztrials(iSplit + (iMod-1)*params.nSplits,:),'LineWidth',0.5);
                        %                 line([vecTimes(:)';vecTimes(:)'],[intTrial*ones(1,numel(vecTimes))-0.5;intTrial*ones(1,numel(vecTimes))+0.5],'Color',params.colors_ztrials{i},'LineWidth',2);
                    end
                    plot([0 0],[0 length(trialselec)],'k--','LineWidth',2)
                    
                    %set fig props
                    ylim([0 length(trialselec)]);
                    xlim([params.t_pre 1e6]);
                    if iSplit==params.nSplits
                        xlabel('Time from stimulus change (s)','FontSize',12);
                        set(gca,'XTick',[-0.5 0 0.5 1]*1e6,'XTickLabel',[-0.5 0 0.5 1],'YTick',[],'FontSize',12)
                    else
                        xlabel('')
                        set(gca,'box','off','XTickLabel',[],'XTick',[],'YTickLabel',[],'YTick',[])
                    end
                    ylabel(params.labels_ztrials{iSplit + (iMod-1)*params.nSplits},'FontSize',12);
                    
                end
            end
            
            tightfig();
            if isfield(params,'exportfig') && params.exportfig
                export_fig(fullfile(params.savedir,sprintf('ExNeuron_Raster_%s',spikeData.cell_ID{cell_idx})),'-eps','-nocrop')
%                 export_fig(fullfile(params.savedir,sprintf('ExNeuron_Raster_%s',spikeData.cell_ID{cell_idx})),'-png','-nocrop')
            end
            
        end
    end
    
end


end