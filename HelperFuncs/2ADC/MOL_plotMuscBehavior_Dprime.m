function MOL_plotMuscBehavior_Dprime(params,dVis_sal,dAud_sal,dVis_musc,dAud_musc,Mice_Sal,Mice_Musc)
%Dprime figure:
% Make figure:
figure; set(gcf,'color','w','units','normalized','Position', [rand()*0.7+0.05 0.4 .15 .45]); hold all;

xpos = [1 2];
fprintf('Visual thr change, Sal vs Musc:\n')
MOL_subPlot_Dprime_Condition(xpos,params.colors_visual_opto{2},dVis_sal(:,1),dVis_musc(:,1),Mice_Sal,Mice_Musc);
% xpos = [3 4];
fprintf('Visual max change UST, Sal vs Musc:\n')
% MOL_subPlot_Dprime_Condition(xpos,params.colors_visual_opto{2},squeeze(dVis_sal(:,1,[1 3])));
MOL_subPlot_Dprime_Condition(xpos,params.colors_visual_opto{1},dVis_sal(:,2),dVis_musc(:,2),Mice_Sal,Mice_Musc);

%Auditory part;
xpos = [3 4];
fprintf('Auditory thr change, Sal vs Musc:\n')
MOL_subPlot_Dprime_Condition(xpos,params.colors_audio_opto{2},dAud_sal(:,1),dAud_musc(:,1),Mice_Sal,Mice_Musc);
% xpos = [7 8];
fprintf('Auditory max change UST, Sal vs Musc:\n')
% MOL_subPlot_Dprime_Condition(xpos,params.colors_visual_opto{2},squeeze(dVis_sal(:,1,[1 3])));
MOL_subPlot_Dprime_Condition(xpos,params.colors_audio_opto{1},dAud_sal(:,2),dAud_musc(:,2),Mice_Sal,Mice_Musc);

%Make up:
ylabel('Dprime')
XTickLabels = repmat({'Sal' 'Mus'},1,2);
set(gca,'XTick',1:length(XTickLabels),'XTickLabels',XTickLabels,'XTickLabelRotation',60);
set(gca,'YTick',[0 1 2])
ylim([0 2])
xlim([0.4 length(XTickLabels)+0.2])
grid on 


end


function MOL_subPlot_Dprime_Condition(xpos,color,dMat1,dMat2,Mice_Sal,Mice_Musc)

y_mean          = [nanmean(dMat1) nanmean(dMat2)];
y_error(1)      = nanstd(dMat1) / sqrt(sum(~isnan(dMat1)));
y_error(2)      = nanstd(dMat2) / sqrt(sum(~isnan(dMat2)));

errorbar(xpos, y_mean,y_error,'-','Color','k','MarkerSize',0.001,'MarkerEdgeColor','k','MarkerFaceColor','k','LineWidth',3);
errorbar(xpos, y_mean,[0 0],[0 0],'o','Color',color,'MarkerSize',12,'MarkerEdgeColor',color,'MarkerFaceColor',color,'LineWidth',0.0001);

%Statistical testing:    signed rank test with bonferroni correction:
tbl             = table([dMat1; dMat2],[ones(size(Mice_Sal)); ones(size(Mice_Sal))*2],[Mice_Sal; Mice_Musc],'VariableNames',{'Hitrate','Treatment','Mouse'}); %Create table for mixed model
lme             = fitlme(tbl,'Hitrate~Treatment+(1|Mouse)'); %construct linear mixed effects model with fixed effect of treatment and random intercept for different mice
stats           = dataset2table(anova(lme,'DFMethod','Satterthwaite')); %Perform ANOVA on model and output as matrix
fprintf('F(%d,%2.0f)=%1.2f, p=%1.3f; \n',stats{2,3},stats{2,4},stats{2,2},stats{2,5})
p = stats{2,5};
if p<0.05
    sigstar(xpos, p) %use sigstar function to identify significance
end
    
% [p] = ranksum(dMat1,dMat2); 
% % p = min([p*3 1]); %bonferroni correction
% if p<0.05
%     sigstar(xpos, p) %use sigstar function to identify significance 
% end

end
