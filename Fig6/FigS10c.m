%% Oude Lohuis et al. 2023 Nat Neurosci
% Triple dissociation of auditory, visual, and motor processing in primary visual cortex
% MOL (C) 2023

% Plot raster of licks and rewards for one sessions
% MOL_SessionRast

params.savedir              = 'E:\Documents\PhD\Figures\Project CHDET\Results - auV1\8ConflictBehavior';
params.exSession            = '2018-01-27_16-18-00';

%% 
[Data] = MOL_GetData('E:','CHDET',{'BehaviorConflict'},{},{params.exSession},{'sessionData' 'trialData'});
sessionData     = Data.sessionData;
trialData       = Data.trialData;

nTrials = length(trialData.session_ID);
for iTrial = 1:nTrials-1
    trialData.lickTime{iTrial} = [trialData.lickTime{iTrial:iTrial+1}];
    trialData.lickSide{iTrial} = [trialData.lickSide{iTrial:iTrial+1}];
end

%% Make figure:
set(0,'defaultAxesFontSize',10)

states          = {'iti'            'stim'       };
statecolors     = {[0.95 0.95 0.95],[0.1 0.7 0.3]};

trialtypes          = {'P'          'X'      'Y'       'C'};
trialcolors         = {[0.5 0.5 0.5],[0.5 0.5 0.95],[0.95 0.5 0.5],[0.7 0.14 0.7]};
params.triallabels  = {'Catch'   'Visual' 'Auditory' 'Conflict'};

if sessionData.VisualLeftCorrectSide
    LeftRewColor       = [0.15 0 0.58];
    RightRewColor      = [0.83 0 0.16];
else
    LeftRewColor       = [0.83 0 0.16];
    RightRewColor      = [0.15 0 0.58];
end


AudioLickColor      = [1 0.47 0.37];
VisualLickColor     = [0 0.67 0.94];

SortBy              = 'stimCange'; %Sort the trials to get a clean overview
AlignOn             = 'stimChange'; %On which timestamp to align as t=0

nTypes              = length(trialtypes);

params.facealpha = 0.5;
h = [];
for sesid = unique(sessionData.session_ID)'
    
    figure;
    set(gcf,'units','normalized','Position',[0.2 0.37 0.52 0.25],'color','w')
    for iType = 1:nTypes
        [temptrialData] = MOL_getTempPerSes(sesid,trialData);
        idx = strcmp(temptrialData.trialType,trialtypes{iType});
        trialFields         = fieldnames(temptrialData);
        for iF = 1:length(trialFields)
            temptrialData.(trialFields{iF}) = temptrialData.(trialFields{iF})(idx);
        end
        
        h(iType) = subplot(1,nTypes,iType);
        
        %Sort the trials
        if isfield(temptrialData,SortBy)
            Delay = temptrialData.(SortBy) - temptrialData.trialStart;
            [~,trialorder] = sort(Delay);
        else
            [~,trialorder] = sort(temptrialData.trialNum);
        end
        
        p = patch([0 1.5 1.5 0],[1 1 length(trialorder)+1 length(trialorder)+1],trialcolors{iType});
        set(p,'FaceAlpha',params.facealpha,'EdgeColor','none');
        
        for i = 1:length(trialorder)
            trial = trialorder(i);
            
            for state = 1:length(states)
                if isfield(temptrialData,strcat(states{state},'Start')) && isfield(temptrialData,strcat(states{state},'End'))
                    
                    ts_start    = temptrialData.(strcat(states{state},'Start'))(trial) - temptrialData.(AlignOn)(trial);
                    ts_end      = temptrialData.(strcat(states{state},'End'))(trial)  - temptrialData.(AlignOn)(trial);
                    
                    if ~isnan(ts_start) && ~isnan(ts_end)
                        if strcmp(states{state},'stim')
                            idx     = strcmp(trialtypes,temptrialData.trialType{trial});
                            ts_end      = 2e6;
                            rectangle('Position',[ts_start/1e6 i (ts_end-ts_start)/1e6 1],'FaceColor',trialcolors{idx},'EdgeColor','none'); hold all;
                        else
                            rectangle('Position',[ts_start/1e6 i (ts_end-ts_start)/1e6 1],'FaceColor',statecolors{state},'EdgeColor','none'); hold all;
                        end
                    end
                end
            end
            
            if isfield(temptrialData,'rewardTime')
                for iRew = 1:length(temptrialData.rewardTime{trial})
                    switch temptrialData.rewardSide{trial}(iRew)
                        case 'L'
                            rectangle('Position',[temptrialData.rewardTime{trial}(iRew)/1e6-temptrialData.(AlignOn)(trial)/1e6 i 0.2 1],'FaceColor',LeftRewColor,'EdgeColor',LeftRewColor); hold all;
                        case 'R'
                            rectangle('Position',[temptrialData.rewardTime{trial}(iRew)/1e6-temptrialData.(AlignOn)(trial)/1e6 i 0.2 1],'FaceColor',RightRewColor,'EdgeColor',RightRewColor); hold all;
                    end
                end
            end
           
            licktimes = [temptrialData.lickTime{trial}];
            licksides = [temptrialData.lickSide{trial}];
            for lick = 1:length(licktimes)
                if sessionData.VisualLeftCorrectSide
                    RightLickColor = AudioLickColor; LeftLickColor = VisualLickColor;
                else RightLickColor = VisualLickColor; LeftLickColor = AudioLickColor;
                end
                licktime = licktimes(lick) - temptrialData.(AlignOn)(trial);
                if ~isnan(licktime) %Sometimes last trial of a session did not get TTL of stimChange
                    if licksides(lick) == 'R'
                        rectangle('Position',[licktime/1e6 i 0.005 1],'FaceColor',RightLickColor,'EdgeColor',RightLickColor,'LineWidth',2); hold all;
                    elseif licksides(lick) == 'L'
                        rectangle('Position',[licktime/1e6 i 0.005 1],'FaceColor',LeftLickColor,'EdgeColor',LeftLickColor,'LineWidth',2); hold all;
                    end
                end
            end
        end
        
        plot([0 0],[1 i+1],'k','LineWidth',0.5);
        
        xmax = max(temptrialData.trialEnd - temptrialData.trialStart);
        if xmax>15; xmax = 15; end
        xlim([-8 1.5])
        ylim([1 i+1])
%         ylabel('Trials')
        title(params.triallabels{iType})
        xlabel('Time relative to stimulus change (s)','FontSize',8)
        set(gca,'Xtick',[-8:2:0 1.5])
        set(gca,'YTick',length(temptrialData.session_ID),'FontSize',10);
        %     set(gca,'Xtick',2:2:floor(xmax),'Ytick',trialorder'+0.5,'yticklabel',cellstr(num2str(trialorder'))')
    end
end
tightfig();

%% For diff trial types: 
for iType = 1:nTypes
    subplot(h(iType));
    F=getframe;
    cla reset
    imagesc(F.cdata);
    set(gca,'XColor','none','YColor','none')
end
filename = sprintf('Conflict_SessionRast_%s.eps',sessionData.session_ID{1});
export_fig(fullfile(params.savedir,filename),gcf)

%% ha
filename = sprintf('Conflict_SessionRast_%s.bmp',sessionData.session_ID{1});
export_fig(fullfile(params.savedir,filename),gcf)
hf = get(h, 'children');
delete(hf{1}); delete(hf{2}); delete(hf{3}); delete(hf{4});

filename = sprintf('Conflict_SessionRast_%s.eps',sessionData.session_ID{1});
export_fig(fullfile(params.savedir,filename),gcf)

%% Legend only:
figure;
rectangle('Position',[0 1 0.2 1],'FaceColor',LeftRewColor,'EdgeColor',LeftRewColor); hold all;
rectangle('Position',[0 3 0.2 1],'FaceColor',RightRewColor,'EdgeColor',RightRewColor); hold all;
rectangle('Position',[0 5 0.005 1],'FaceColor',RightLickColor,'EdgeColor',RightLickColor,'LineWidth',0.8); hold all;
rectangle('Position',[0 7 0.005 1],'FaceColor',LeftLickColor,'EdgeColor',LeftLickColor,'LineWidth',0.8); hold all;
xlim([-.4 .4])
ylim([0 9])
filename = sprintf('Conflict_SessionRast_Legend.eps');
export_fig(fullfile(params.savedir,filename),gcf)
