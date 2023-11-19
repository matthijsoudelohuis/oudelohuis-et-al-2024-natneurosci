function MOL_plotTuning_OriFreqTuning(sessionData,trialData,spikeData,cell_IDs,params)

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
            
            %Compute histogram:
            events_ts               = temptrialData.(params.AlignOn);
            spikes_ts               = spikeData.ts{cell_idx};
            hist_mat                = calc_psth(events_ts,spikes_ts,params);    %Construct histogram matrix:
            resp_mat                = calc_resp_from_psth(hist_mat,params);
            
            figure; set(gcf,'units','normalized','Position',[0.1 0.5 0.2 0.17],'color','w'); hold all;
            
            for iMod = 1:2
                meantoplot = [];
                for iSplit = 1:params.nSplits %Store the mean response for each of these splits
                    meantoplot(iSplit) = nanmean(resp_mat(splits{iSplit,iMod}),1); %#ok<*AGROW>
                    errortoplot(iSplit) = nanstd(resp_mat(splits{iSplit,iMod}),1) / sqrt(sum(splits{iSplit,iMod})); %#ok<*AGROW>
                end
                
                subplot(1,2,iMod); hold all;
                errorbar(meantoplot,errortoplot,'k','MarkerSize',0.1,'LineWidth',0.5)
                for iSplit = 1:params.nSplits 
                    errorbar(iSplit,meantoplot(iSplit),0,'.','Color',params.colors_ztrials(iSplit + (iMod-1)*params.nSplits,:),'LineWidth',0.1,'MarkerSize',25);
                end
                
                if iMod==1
                    title('Orientation','FontSize',12)
                elseif iMod==2
                    title('Frequency','FontSize',12)
                end
                
                xlabel('')
                xlim([0.5 8.5]);
                set(gca,'XTickLabel',params.labels_ztrials((1:8) + (iMod-1)*params.nSplits))
                set(gca,'XTick',1:8,'XTickLabelRotation',45,'TickDir','out')
                ylabel('')
                ylim_sub(iMod,:) = get(gca,'YLim');
            end
            
            y_max = round(max(max(ylim_sub)));
            for iMod = 1:2
                subplot(1,2,iMod);
                set(gca,'YLim',[0 y_max],'YTick',[0 y_max])
            end

            
%             tightfig();
            if isfield(params,'exportfig') && params.exportfig
                export_fig(fullfile(params.savedir,sprintf('ExNeuron_Tuning_%s',spikeData.cell_ID{cell_idx})),'-eps','-nocrop')
%                 export_fig(fullfile(params.savedir,sprintf('ExNeuron_Tuning_%s',spikeData.cell_ID{cell_idx})),'-png','-nocrop')
            end
            
        end
    end
    
end


end