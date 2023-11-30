%% Load exported orientation decoding scores from python: 
datadir         = 'E:\Matlab\oudelohuis-et-al-2024-natneurosci\FigS3\';
filename        = fullfile(datadir,'SourceData_FigS3_audiofreq_in_big_change_feb8_singlesess_200.csv');
tbl             = readtable(filename);
writetable(tbl,fullfile(datadir,'SourceData_FigS3_audiofreq_in_big_change.xlsx'));

%% Correlation frequency decoding performance spikes and video:
lme             = fitlme(tbl,'score_video~score_spikes+(1|animal_id)'); %construct linear mixed effects model with fixed effect of temporal window and random intercept for different mice
stats           = dataset2table(anova(lme,'DFMethod','Satterthwaite')); %Perform ANOVA on model and output as matrix
fprintf('Correlation frequency decoding performance spikes and video:\n')
fprintf('R=%1.2f,F(%d,%2.0f)=%1.2f, p=%1.3e; \n',sqrt(lme.Rsquared.Ordinary),stats{2,3},stats{2,4},stats{2,2},stats{2,5})

%% Difference frequency decoding performance spikes and video:

Y_score     = [tbl.score_spikes; tbl.score_video];
X_met       = [ones(size(tbl,1),1); ones(size(tbl,1),1)*2];
G_mou       = [tbl.animal_id; tbl.animal_id];
tbl2            = table(Y_score,X_met,G_mou,'VariableNames',{'score','method','animal_id'}); %Create table for mixed model

lme             = fitlme(tbl2,'score~method+(1|animal_id)'); %construct linear mixed effects model with fixed effect of temporal window and random intercept for different mice
stats           = dataset2table(anova(lme,'DFMethod','Satterthwaite')); %Perform ANOVA on model and output as matrix
fprintf('Difference freq decoding performance spikes and video:\n')
fprintf('F(%d,%2.0f)=%1.2f, p=%1.3f; \n',stats{2,3},stats{2,4},stats{2,2},stats{2,5})

%% Load exported orientation decoding scores from python: 
filename        = fullfile(datadir,'SourceData_FigS3_visualori_in_big_change_feb8_singlesess_200.csv');
tbl             = readtable(filename);
writetable(tbl,fullfile(datadir,'SourceData_FigS3_visualori_in_big_change.xlsx'));

%% Correlation orientation decoding performance spikes and video::
lme             = fitlme(tbl,'score_video~score_spikes+(1|animal_id)'); %construct linear mixed effects model with fixed effect of temporal window and random intercept for different mice
stats           = dataset2table(anova(lme,'DFMethod','Satterthwaite')); %Perform ANOVA on model and output as matrix
fprintf('Correlation orientation decoding performance spikes and video:\n')
fprintf('R=%1.2f,F(%d,%2.0f)=%1.2f, p=%1.2f; \n',sqrt(lme.Rsquared.Ordinary),stats{2,3},stats{2,4},stats{2,2},stats{2,5})
p = stats{2,5};

%% Difference orientation decoding performance spikes and video:
Y_score         = [tbl.score_spikes; tbl.score_video];
X_met           = [ones(size(tbl,1),1); ones(size(tbl,1),1)*2];
G_mou           = [tbl.animal_id; tbl.animal_id];
tbl2            = table(Y_score,X_met,G_mou,'VariableNames',{'score','method','animal_id'}); %Create table for mixed model

lme             = fitlme(tbl2,'score~method+(1|animal_id)'); %construct linear mixed effects model with fixed effect of temporal window and random intercept for different mice
stats           = dataset2table(anova(lme,'DFMethod','Satterthwaite')); %Perform ANOVA on model and output as matrix
fprintf('Difference ori decoding performance spikes and video:\n')
fprintf('F(%d,%2.0f)=%1.2f, p=%1.3e; \n',stats{2,3},stats{2,4},stats{2,2},stats{2,5})

