function plotRaster_AuResp(sessionData,trialData,spikeData,cell_IDs,params)

figure; set(gcf,'units','normalized','Position',[0.1 0.4 0.5 0.4],'color','w'); hold all;

nNeurons                = length(cell_IDs);

params.nSplits          = 1;

params.AlignOn          = 'stimChange';

for iNeuron = 1:nNeurons %Loop over all neurons:
    
    cell_idx = strcmp(spikeData.cell_ID,cell_IDs(iNeuron));
    
    if any(cell_idx)
        %Get the relevant data for each session individually:
        temptrialData        = MOL_getTempPerSes(spikeData.session_ID(cell_idx),trialData);
        
        temptrialData.visualOriPostChangeNorm = round(temptrialData.visualOriPostChangeNorm/2);
        temptrialData.audioFreqPostChangeNorm = round(temptrialData.audioFreqPostChangeNorm/2);
        
        splits          = {};
%         splits{1}       = ismember(temptrialData.trialType,{'Y'}) & temptrialData.audioFreqChangeNorm==3;
        splits{1}       = temptrialData.audioFreqChangeNorm==3;
        %         splits{1}       = ismember(temptrialData.trialType,{'Y'}) & temptrialData.audioFreqChangeNorm==2;
        %         splits{2}       = ismember(temptrialData.trialType,{'X'}) & temptrialData.visualOriChangeNorm==3 & temptrialData.visualOriPostChangeNorm==2;
        %         splits{3}       = ismember(temptrialData.trialType,{'Y'}) & temptrialData.audioFreqChangeNorm==3 & temptrialData.audioFreqPostChangeNorm==1;
        %         splits{4}       = ismember(temptrialData.trialType,{'Y'}) & temptrialData.audioFreqChangeNorm==3 & temptrialData.audioFreqPostChangeNorm==2;
        
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
        
        timeidx                 = params.xtime>0 & params.xtime<0.3e6;
        if  nanmean(nanmean(hist_mat(splits{1} & temptrialData.audioFreqPostChangeNorm==1,timeidx),1),2) > nanmean(nanmean(hist_mat(splits{1} & temptrialData.audioFreqPostChangeNorm==2,timeidx),1),2)
            splits{1}       = splits{1} & temptrialData.audioFreqPostChangeNorm==1;
        else splits{1}       = splits{1} & temptrialData.audioFreqPostChangeNorm==2;
        end
        
        for iSplit = 1:params.nSplits %Store the mean response for each of these splits
            ratemat(iSplit,:) = nanmean(hist_mat(splits{iSplit},:),1); %#ok<*AGROW>
            
%             [h,templat,tempz,tempr]     = getZeta(spikes_ts*1e-6,events_ts(splits{iSplit})*1e-6 - 0.3,1,[],0,4,[0 0.2],[]);
%             tempr.vecT = (tempr.vecT - 0.3) * 1e6;
            
        end
        
        subplot(2,nNeurons,iNeuron); hold all;
        for iSplit = 1:params.nSplits %Store the mean response for each of these splits
%                         handles(iSplit) = plot(params.xtime,ratemat(iSplit,:),'-','Color',params.colors_ztrials{iSplit},'LineWidth',3);
                        handles(iSplit) = plot(params.xtime,ratemat(iSplit,:),'-','Color','k','LineWidth',1.5);
%             handles(iSplit) = plot(tempr.vecT,tempr.vecRate,'-','Color','k','LineWidth',1.5);
        end
        
        xlabel('')
        xlim([params.t_pre params.t_post]);
        ylim([0 ceil(max(max(ratemat)))])
        set(gca,'XTick',[],'YTick',get(gca,'YLim'))
        plot([0 0],[0 max(max(ratemat))],'--','Color',[0.6 0.6 0.6],'LineWidth',1)

        if iNeuron==1
            ylabel('Firing rate (sp s-1)','FontSize', 15)
        else ylabel('')
        end
        
        for i = 1:params.nSplits
            subplot(2,nNeurons,iNeuron+nNeurons); hold all;
            
            %Sort trials based on post-change orientation and response latency (orientation is sorted by adding a lot of latency)
            trialselec          = find(splits{i}');
            temp                = temptrialData.responseLatency(trialselec);
            temp(isnan(temp))   = randn(sum(isnan(temp)),1);
            [~,idx]             = sort(temp + 1e13*temptrialData.audioFreqPostChangeNorm(trialselec),'descend');
            trialselec          = trialselec(idx);
            
            for intTrial = 1:length(trialselec)%find(trialselec') %1:numel(vecTrialStarts)
                
                vecTimes = vecTimePerSpike(vecTrialPerSpike==trialselec(intTrial));
                lat = temptrialData.responseLatency(trialselec(intTrial));
                if ~isnan(lat)
                    line([lat lat],[intTrial-0.5;intTrial+0.5],'Color','k','LineWidth',2);
                end
                line([vecTimes(:)';vecTimes(:)'],[intTrial*ones(1,numel(vecTimes))-1;intTrial*ones(1,numel(vecTimes))+1],'Color',params.colors_experiments{iNeuron},'LineWidth',1);
%                 line([vecTimes(:)';vecTimes(:)'],[intTrial*ones(1,numel(vecTimes))-1;intTrial*ones(1,numel(vecTimes))+1],'Color','k','LineWidth',1);
                %                 line([vecTimes(:)';vecTimes(:)'],[intTrial*ones(1,numel(vecTimes))-0.5;intTrial*ones(1,numel(vecTimes))+0.5],'Color',params.colors_ztrials{i},'LineWidth',2);
            end
            plot([0 0],[0 length(trialselec)],'--','Color',[0.6 0.6 0.6],'LineWidth',1)
            
            %set fig props
            ylim([0 length(trialselec)]);
            xlim([params.t_pre params.t_post]);
            set(gca,'XTick',[-0.25 0 0.25 0.5]*1e6,'XTickLabel',[-250 0 250 500],'YTick',[])
            set(gca,'box','off')
            if iNeuron==2
                xlabel('Time from stimulus change (s)');
            else
                xlabel('')
            end
            ylabel(params.labels_neurons{iNeuron},'FontSize',15);
        end
    end
end

tightfig();
if isfield(params,'exportfig') && params.exportfig
    export_fig(fullfile(params.savedir,'PSTH_ExNeurons_3cohorts'),'-eps','-nocrop')
    %                 export_fig(fullfile(params.savedir,sprintf('ExNeuron_%s',spikeData.cell_ID{cell_idx})),'-png','-nocrop')
end

end