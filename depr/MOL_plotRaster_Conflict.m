function MOL_plotRaster_Conflict(sessionData,trialData,spikeData,cell_IDs,params)

nNeurons                = length(cell_IDs);

params.trialcategories = 'Confl_Stim';
params.nSplits          = 8;
params.labels_ztrials   = {'' 'a-' 'A-'; '-v' 'av' 'Av'; '-V' 'aV' 'AV'};
postclass = 2;

%         params.trialcategories = 'Confl_Choice';
% params.nSplits          = 4;

params.AlignOn          = 'stimChange';

for iNeuron = 1:nNeurons %Loop over all neurons:
    
    cell_idx = strcmp(spikeData.cell_ID,cell_IDs(iNeuron));
    
    if any(cell_idx)
        %Get the relevant data for each session individually:
        temptrialData        = MOL_getTempPerSes(spikeData.session_ID(cell_idx),trialData);
        
        switch params.trialcategories
            case 'Confl_Stim'
                idx_C           = strcmp(temptrialData.trialType,'C');
                                
                if numel(unique(temptrialData.visualOriPostChangeNorm(idx_C)))==1
                    splits          = {};
                    splits{1,1}     = [];
                    splits{2,1}     = ismember(temptrialData.trialType,{'X'}) & temptrialData.visualOriChangeNorm==2 & temptrialData.visualOriPostChangeNorm==unique(temptrialData.visualOriPostChangeNorm(idx_C));
                    splits{3,1}     = ismember(temptrialData.trialType,{'X'}) & temptrialData.visualOriChangeNorm==3 & temptrialData.visualOriPostChangeNorm==unique(temptrialData.visualOriPostChangeNorm(idx_C));
                    splits{1,2}     = ismember(temptrialData.trialType,{'Y'}) & temptrialData.audioFreqChangeNorm==2 & temptrialData.audioFreqPostChangeNorm==unique(temptrialData.audioFreqPostChangeNorm(idx_C));
                    splits{1,3}     = ismember(temptrialData.trialType,{'Y'}) & temptrialData.audioFreqChangeNorm==3 & temptrialData.audioFreqPostChangeNorm==unique(temptrialData.audioFreqPostChangeNorm(idx_C));
                    
                    splits{2,2}     = ismember(temptrialData.trialType,{'C'}) & temptrialData.visualOriChangeNorm==2 & temptrialData.audioFreqChangeNorm==2;
                    splits{2,3}     = ismember(temptrialData.trialType,{'C'}) & temptrialData.visualOriChangeNorm==2 & temptrialData.audioFreqChangeNorm==3;
                    splits{3,2}     = ismember(temptrialData.trialType,{'C'}) & temptrialData.visualOriChangeNorm==3 & temptrialData.audioFreqChangeNorm==2;
                    splits{3,3}     = ismember(temptrialData.trialType,{'C'}) & temptrialData.visualOriChangeNorm==3 & temptrialData.audioFreqChangeNorm==3;
                else
                    temptrialData.visualOriPostChangeNorm = round(temptrialData.visualOriPostChangeNorm/2);
                    temptrialData.audioFreqPostChangeNorm = round(temptrialData.audioFreqPostChangeNorm/2);
                    
                    splits          = {};
                    splits{1,1}     = [];
                    splits{2,1}     = ismember(temptrialData.trialType,{'X'}) & temptrialData.visualOriChangeNorm==2 & temptrialData.visualOriPostChangeNorm==postclass;
                    splits{3,1}     = ismember(temptrialData.trialType,{'X'}) & temptrialData.visualOriChangeNorm==3 & temptrialData.visualOriPostChangeNorm==postclass;
                    splits{1,2}     = ismember(temptrialData.trialType,{'Y'}) & temptrialData.audioFreqChangeNorm==2 & temptrialData.audioFreqPostChangeNorm==postclass;
                    splits{1,3}     = ismember(temptrialData.trialType,{'Y'}) & temptrialData.audioFreqChangeNorm==3 & temptrialData.audioFreqPostChangeNorm==postclass;
                    
                    splits{2,2}     = ismember(temptrialData.trialType,{'C'}) & temptrialData.visualOriChangeNorm==2 & temptrialData.audioFreqChangeNorm==2 & temptrialData.visualOriPostChangeNorm==postclass & temptrialData.audioFreqPostChangeNorm==postclass;
                    splits{2,3}     = ismember(temptrialData.trialType,{'C'}) & temptrialData.visualOriChangeNorm==2 & temptrialData.audioFreqChangeNorm==3 & temptrialData.visualOriPostChangeNorm==postclass & temptrialData.audioFreqPostChangeNorm==postclass;
                    splits{3,2}     = ismember(temptrialData.trialType,{'C'}) & temptrialData.visualOriChangeNorm==3 & temptrialData.audioFreqChangeNorm==2 & temptrialData.visualOriPostChangeNorm==postclass & temptrialData.audioFreqPostChangeNorm==postclass;
                    splits{3,3}     = ismember(temptrialData.trialType,{'C'}) & temptrialData.visualOriChangeNorm==3 & temptrialData.audioFreqChangeNorm==3 & temptrialData.visualOriPostChangeNorm==postclass & temptrialData.audioFreqPostChangeNorm==postclass;
                end
                
                
            case 'Confl_Choice'
                splits          = {};
                splits{1}       = ismember(temptrialData.trialType,{'X'}) & temptrialData.vecResponse==2;
                splits{2}       = ismember(temptrialData.trialType,{'Y'}) & temptrialData.vecResponse==1;
                
                splits{3}       = ismember(temptrialData.trialType,{'C'}) & temptrialData.vecResponse==2;
                splits{4}       = ismember(temptrialData.trialType,{'C'}) & temptrialData.vecResponse==1;
        end

        if all(cellfun(@sum,splits([0 1 1; 1 1 1; 1 1 1]==1))>=params.minNtrials)
            
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
            
            ratelim = [1000 0]; %define lower and higher limits by going through all combinations and finding max and min and updating lim vector
            for i = 1:3
                for j = 1:3
                    ratelim(1) = min([ratelim(1) min(nanmean(hist_mat(splits{i,j},:),1))*0.9]);
                    ratelim(2) = max([ratelim(2) max(nanmean(hist_mat(splits{i,j},:),1))*1.1]);
                end
            end
            ratelim = round(ratelim);
            
            figure; set(gcf,'units','normalized','Position',[0.1 0.4 0.26 0.4],'color','w'); hold all;
            
            for i = 1:3 %for each visual condition
                for j = 1:3 %for each auditory condition
                    if ~isempty(splits{i,j}) %if there is a selection of trials (not for 1,1)
                        subplot(3,3,j+(i-1)*3,'position',[0.05+(j-1)*0.3 0.95-i*0.3 0.27 0.27]); hold all;
                        set(gca, 'color', [0.96 0.96 0.96])
                        
                        %Sort trials based on response latency                      
                        trialselec          = find(splits{i,j}');
                        temp                = temptrialData.responseLatency(trialselec);
                        temp(isnan(temp))   = 0;
                        [~,idx]             = sort(temp,'descend');
                        trialselec          = trialselec(idx);
                        
                        for intTrial = 1:length(trialselec) %start plotting ticks for each trial
                            
                            vecTimes = vecTimePerSpike(vecTrialPerSpike==trialselec(intTrial));
                            lat = temptrialData.responseLatency(trialselec(intTrial));
                            if ~isnan(lat)
                                line([lat lat],[intTrial-1;intTrial],'Color','k','LineWidth',1);
                            end
                            line([vecTimes(:)';vecTimes(:)'],[intTrial*ones(1,numel(vecTimes))-1;intTrial*ones(1,numel(vecTimes))],'Color',params.colors_ticks{i,j},'LineWidth',0.15);
                            %                 line([vecTimes(:)';vecTimes(:)'],[intTrial*ones(1,numel(vecTimes))-0.5;intTrial*ones(1,numel(vecTimes))+0.5],'Color',params.colors_ticks{i},'LineWidth',2);
                        end
                        plot([0 0],[0 length(trialselec)],'k--','LineWidth',2)
                        
                        %set fig props (make up)
                        ylim([0 length(trialselec)]);
                        xlim([params.t_pre params.t_post]);
                        if i==3 && j==2
                            xlabel('Time from stimulus change (s)');
                            set(gca,'box','off','XTickLabel',[],'XTick',[-0.5 0 0.5 1]*1e6,'XTickLabel',[-0.5 0 0.5 1],'YTick',[])
                        else
                            xlabel('')
                            set(gca,'box','off','XTickLabel',[],'XTick',[],'YTickLabel',[],'YTick',[])
                        end
                        title(params.labels_ztrials{i,j},'FontSize',15); %make title based on av combination
                        
                    end
                    
                    %Plot the firing rate through the figure
                    ratemat                 = nanmean(hist_mat(splits{i,j},:),1); %#ok<*AGROW>
                    ratemat                 = ratemat - ratelim(1);
                    ratemat                 = ratemat / ratelim(2) * sum(splits{i,j});
                    plot(params.xtime,ratemat,'-','Color',[0.2 0.2 0.2],'LineWidth',1.5);
                    if i==3 && j==3
                        yyaxis('right');
                        ylim(ratelim)
                        set(gca,'YTick',ratelim)
                    end

                end
            end
            
            h = subplot(3,3,1);
            set(gca,'XTick',[],'YTick',[],'box','off')
            h.Color = 'none';
            title(cell_IDs(iNeuron),'FontSize',12)
            
%             tightfig();
            if isfield(params,'exportfig') && params.exportfig
                export_fig(fullfile(params.savedir,sprintf('ExNeuron_%s',spikeData.cell_ID{cell_idx})),'-eps','-nocrop')
%                 export_fig(fullfile(params.savedir,sprintf('ExNeuron_%s',spikeData.cell_ID{cell_idx})),'-png','-nocrop')
            end
        end
    end
    
end


end
