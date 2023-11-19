function MOL_plotMuscBehavior_Rates(params,TotalResp_Sal,TotalResp_Musc,Mice_Sal,Mice_Musc)

figure; set(gcf,'color','w','units','normalized','Position', [0.4 0.38 0.25 0.27]);    hold all;

% Settings:
params.auprobepos       = 0.001;
params.auticks          = [1/256 1/64 1/8 1/2];
params.auticklabels     = {'Probe' '1/256' '1/64' '1/8' '1/2'};
params.auxaxislabel     = 'Change in octave';
params.auystatslabel    = 'Auditory threshold (partial octave)';

params.visprobepos      = 0.5;
params.visticks         = [2 5 15 30 90];
params.vistickslabels   = ['Probe' num2cell(params.visticks)];
params.visxaxislabel    = 'Change in orientation (Degrees)';
params.visystatslabel   = 'Visual threshold (Degrees)';

params.yticks           = [0 0.25 0.5 0.75 1];

% Construct x from position of the probe trial and conditions
xdata_vis               = [params.visprobepos 12 90];

% Construct x from position of the probe trial and conditions
xdata_au                = [params.auprobepos 1/32 1/2];

params.auticks         = [1/32 1/2];
params.auticklabels   = {'Catch' 'Thr' 'Max'};

params.visticks         = [12 90];
params.vistickslabels   = {'Catch' 'Thr' 'Max'};

%Audio:
subplot(1,2,1); hold all;
%Get the hit rates for audio conditions in this condition:
%=full first dimension, vis=1 (no change), resp=1,all sessions
datatoplot_sal = squeeze(TotalResp_Sal(:,1,1,:));
params.linehandles(1) = errorbar(xdata_au,nanmean(datatoplot_sal,2),nanstd(datatoplot_sal,[],2)/sqrt(size(TotalResp_Sal,4)),'.-',...
    'Color',params.colors_audio_opto{1},'LineWidth',1,'MarkerSize',35);

datatoplot_mus = squeeze(TotalResp_Musc(:,1,1,:));
params.linehandles(2) = errorbar(xdata_au,nanmean(datatoplot_mus,2),nanstd(datatoplot_mus,[],2)/sqrt(size(TotalResp_Musc,4)),'.:',...
    'Color',params.colors_audio_opto{2},'LineWidth',1,'MarkerSize',35);

fprintf('Auditory hit rates: \n',numel(datatoplot_sal(1,:)),numel(datatoplot_mus(1,:)))
for i = 1:3
    tbl             = table([datatoplot_sal(i,:)'; datatoplot_mus(i,:)'],[ones(size(Mice_Sal)); ones(size(Mice_Sal))*2],[Mice_Sal; Mice_Musc],'VariableNames',{'Hitrate','Treatment','Mouse'}); %Create table for mixed model
    lme             = fitlme(tbl,'Hitrate~Treatment+(1|Mouse)'); %construct linear mixed effects model with fixed effect of temporal window and random intercept for different mice
    stats           = dataset2table(anova(lme,'DFMethod','Satterthwaite')); %Perform ANOVA on model and output as matrix
    fprintf('%s: F(%d,%2.0f)=%1.2f, p=%1.3f; \n',params.auticklabels{i},stats{2,3},stats{2,4},stats{2,2},stats{2,5})
    p = stats{2,5};
    if p<0.05
        sigstar([xdata_au(i)-0.1 xdata_au(i)+0.1], p) %use sigstar function to identify significance
    end
    writetable(tbl,sprintf('SourceData_Fig5d_AudHitrate_Saliency_%d_Muscimol.xlsx',i))

end

legend(params.linehandles,{'Au - Saline' 'Au - Musc'},'FontSize',10,'Location','NorthEast');
legend boxoff;

%Visual:
subplot(1,2,2); hold all;

%Get the hit rates for visual conditions in this opto condition:
%=au=1 (no change), full visual dimension, iOpto, resp=2,all sessions
datatoplot_sal = squeeze(TotalResp_Sal(1,:,2,:));
if size(datatoplot_sal,1)==1; datatoplot_sal = datatoplot_sal'; end
params.linehandles(1) = errorbar(xdata_vis,nanmean(datatoplot_sal,2),nanstd(datatoplot_sal,[],2)/sqrt(size(TotalResp_Sal,4)),'.-',...
    'Color',params.colors_visual_opto{1},'LineWidth',1,'MarkerSize',35);

datatoplot_mus = squeeze(TotalResp_Musc(1,:,2,:));
if size(datatoplot_mus,1)==1; datatoplot_mus = datatoplot_mus'; end
params.linehandles(2) = errorbar(xdata_vis,nanmean(datatoplot_mus,2),nanstd(datatoplot_mus,[],2)/sqrt(size(TotalResp_Musc,4)),'.:',...
    'Color',params.colors_visual_opto{2},'LineWidth',1,'MarkerSize',35);

fprintf('Visual hit rates: \n')
for i = 1:3
    tbl             = table([datatoplot_sal(i,:)'; datatoplot_mus(i,:)'],[ones(size(Mice_Sal)); ones(size(Mice_Sal))*2],[Mice_Sal; Mice_Musc],'VariableNames',{'Hitrate','Treatment','Mouse'}); %Create table for mixed model
    lme             = fitlme(tbl,'Hitrate~Treatment+(1|Mouse)'); %construct linear mixed effects model with fixed effect of temporal window and random intercept for different mice
    stats           = dataset2table(anova(lme,'DFMethod','Satterthwaite')); %Perform ANOVA on model and output as matrix
    fprintf('%s: F(%d,%2.0f)=%1.2f, p=%1.3f; \n',params.vistickslabels{i},stats{2,3},stats{2,4},stats{2,2},stats{2,5})
    p = stats{2,5};
    if p<0.05
        sigstar([xdata_vis(i)-0.1 xdata_vis(i)+0.1], p) %use sigstar function to identify significance
    end
    
    writetable(tbl,sprintf('SourceData_Fig5d_VisHitrate_Saliency_%d_Muscimol.xlsx',i))

end
fprintf('(F-test LMM model, n=%d,%d sessions)\n ',numel(datatoplot_sal(1,:)),numel(datatoplot_mus(1,:)))

legend(params.linehandles,{'Vis - Saline' 'Vis - Musc'},'FontSize',10,'Location','NorthEast');
legend boxoff;
MOL_Psy2Sided_FigMakeup(params)


end