function MOL_UpdateCSDfig()

global mainfig params temptrialData templfpData

figure(mainfig)

%% Panel 1 and 2: CSD power
% length(splits) = 2;

for iSplit = 1:2
    
    subplot(1,2,iSplit)
    cla;
    [meancsd,meanresp]                  = calcCSD(params,templfpData,temptrialData,iSplit);
    figure(mainfig);
    params.cscale                       = [-max(max(meancsd))*0.9 max(max(meancsd))*0.9];
    
    imagesc(params.xtime,templfpData.sortedChannelDepth,meancsd,params.cscale); hold on;
    xlabel('Time from stimulus (s)','FontSize', 15)
    set(gca,'YDir','normal')
    
    switch params.colormap
        case 'redblue'
            h       = coolwarm + repmat(0.1-abs(linspace(-0.1,.1,64))',1,3);
            h(h>1)  = 1;
            if params.sinkhot
                h = flipud(h);
            end
            colormap(h);
        case 'parula'
            if params.sinkhot
                colormap(flipud(parula));
            else                     colormap(parula);
            end
    end
    
    offsetmat = repmat(templfpData.sortedChannelDepth,1,length(params.xtime));
    plot(params.xtime,meanresp*20e4 + offsetmat,'k','LineWidth',0.5); hold on;
    
    ylabel('Channel','FontSize', 15)
    plot([0 0],ylim,'k','LineWidth',1);
    xlim([params.xtime(1) params.xtime(end)]);
    ylim(params.Ylim_all);
    set(gca,'XTick',params.xticks,'XTickLabels',params.xticklabels)
    xlabel('Time from stimulus (ms)','FontSize', 15)
%     set(gca, 'YTick', flipud(templfpData.sortedChannelDepth), 'YTickLabels', flipud(templfpData.sortedChannelNum),'FontSize', 10)
    set(gca, 'YTick', templfpData.sortedChannelDepth, 'YTickLabels', templfpData.sortedChannelNum,'FontSize', 10)
%     set(gca, 'YTick', flipud(templfpData.sortedChannelDepth), 'YTickLabels', templfpData.sortedChannelNum,'FontSize', 10)
end

end

function [meancsd,meanresp] = calcCSD(params,templfpData,temptrialData,iSplit)

splits          = {};
splits{1}       = ismember(temptrialData.trialType,{'X' 'V'}) & temptrialData.visualOriChangeNorm==3 & ~(temptrialData.hasphotostim==1) & ismember(temptrialData.vecResponse,[1 2]);
splits{2}       = ismember(temptrialData.trialType,{'Y' 'A'}) & temptrialData.audioOctChangeNorm==3 & ~(temptrialData.hasphotostim==1) & ismember(temptrialData.vecResponse,[1 2]);

splits          = {};
%         splits{1}       = ismember(temptrialData.trialType,{'X' 'V'}) & temptrialData.visualOriChangeNorm==3 & temptrialData.correctResponse==0 & ~(temptrialData.hasphotostim==1);
splits{1}       = ismember(temptrialData.trialType,{'X' 'V'}) & temptrialData.visualOriChangeNorm==3 & ~(temptrialData.hasphotostim==1);
%         splits{2}       = ismember(temptrialData.trialType,{'Y' 'A'}) & temptrialData.visualOriChangeNorm==3 & temptrialData.correctResponse==0 & ~(temptrialData.hasphotostim==1);
%         splits{2}       = ismember(temptrialData.trialType,{'Y' 'A'}) & temptrialData.audioFreqChangeNorm==3 & ~(temptrialData.hasphotostim==1);
splits{2}       = ismember(temptrialData.trialType,{'Y' 'A'}) & temptrialData.audioOctChangeNorm==3 & ~(temptrialData.hasphotostim==1);
        
        
%make butter filter:
[params.B_butt,params.A_butt]   = butter(params.ord_butter,params.lp_butter/(templfpData.fs(1)/2));

if params.UseKaiser
    [N,Wn,beta,ftype]   = kaiserord([params.hp_kaiser-params.dp_kaiser params.hp_kaiser params.lp_kaiser params.lp_kaiser+params.dp_kaiser],[0 1 0],[1e-6 1 1e-6],templfpData.fs(1));
    B_kai               = fir1(N,Wn,ftype,kaiser(N+1,beta),'noscale');
    A_kai               = 1;
end

events                  = temptrialData.stimChange(splits{iSplit});

nEvents                         = length(events);
lfpmat                          = zeros(params.nChannels,length(params.xtime),nEvents); %init lfp matrix

%Get lfp per trial and apply low-pass filter:
for ev = 1:nEvents
    
    tempsignal      = templfpData.sortedSignal(:,templfpData.ts>events(ev)+params.t_pre & templfpData.ts<events(ev)+params.t_post);
    lfpmat(:,:,ev)  = tempsignal(:,1:length(params.xtime));
    
    if params.UseButter
        for iCh = 1:params.nChannels
            lfpmat(iCh,:,ev)            = filtfilt(params.B_butt,params.A_butt,lfpmat(iCh,:,ev));
        end
    end
    
    if params.UseKaiser
        for iCh = 1:params.nChannels
            lfpmat(iCh,:,ev)            = filtfilt(B_kai,A_kai,lfpmat(iCh,:,ev));
        end
    end
       
end

% Replace bad channels with mean of adjacent channels
for iCh = find(~templfpData.sortedisgood)'
    chbelow          = find(templfpData.sortedChannelDepth > templfpData.sortedChannelDepth(iCh) & templfpData.sortedisgood,1,'first');
    chabove          = find(templfpData.sortedChannelDepth < templfpData.sortedChannelDepth(iCh) & templfpData.sortedisgood,1,'last');
    lfpmat(iCh,:,:)  = mean(lfpmat([chbelow chabove],:,:),1);
end

%Loop over the different shanks (CSD does not work well with
%considering no horizontal spacing)
for iShank = 1:params.nShanks
    ShankChannelSel                     = params.ChannelSel([1:params.nChannelsShank] + params.nChannelsShank*(iShank-1));
    
    idx                                 = ismember(templfpData.sortedChannelNum,ShankChannelSel);
    sortedShankChannelNum               = templfpData.sortedChannelNum(idx);
    sortedShankChannelDepth             = templfpData.sortedChannelDepth(idx);
    %     sortedShankChannelIdx       = templfpData.sortedChannelNum(ismember(templfpData.sortedChannelNum,ShankChannelSel));
    
    shanklfpmat                         = lfpmat(idx,:,:);

    switch params.csdmethod
        case 'DoubleDiff'
            meancsd                             = csdfromlfp_doublediff(mean(shanklfpmat,3));
        case 'NicholsonFreeman'
            meancsd                             = csdfromlfp_NicholsonFreeman(mean(shanklfpmat,3),templfpData);
    end
   
    %Baseline subtraction of csd:
    meancsd_shank{iShank}               = meancsd - repmat(mean(meancsd(:,params.xtime<-0.025),2),1,length(params.xtime));
    ChannelY_shank{iShank}              = sortedShankChannelDepth;
end

meanresp                = mean(lfpmat(:,:,:),3);
baseline                = repmat(mean(meanresp(:,params.xtime<0),2),1,length(params.xtime));
meanresp                = meanresp - baseline;

meancsd_shank{1}        = tointerpol2(meancsd_shank{1},ChannelY_shank{1}',templfpData.sortedChannelDepth);
meancsd_shank{2}        = tointerpol2(meancsd_shank{2},ChannelY_shank{2}',templfpData.sortedChannelDepth);

meancsd_shank{1}        = meancsd_shank{1}/max(max(meancsd_shank{1}));
meancsd_shank{2}        = meancsd_shank{2}/max(max(meancsd_shank{2}));

meancsd                 = cat(3,meancsd_shank{1},meancsd_shank{2});
meancsd                 = nanmean(meancsd,3);

end

function csd = csdfromlfp_doublediff(lfp)

% Vaknin transform: (Vaknin 1985: add channels above and below assuming isoconductivity):
lfp_vaknin      = [lfp(1,:); lfp(1,:); lfp; lfp(end,:); lfp(end,:)];

% Filter temporally and spatially: is key to this
spat_filter     = fspecial('gaussian',[3 5],1.1);
lfp_vaknin      = conv2(lfp_vaknin,spat_filter,'same');

% compute the CSD map
csd             = diff(lfp_vaknin,2,1); %Taking the double derivative from up to down

% Remove first and last channels of csd map to ensure same size as original
% lfp matrix and set outer channels to zero (nonsensical csd values)
csd             = [zeros(1,size(csd,2));csd(3:end-2,:);zeros(1,size(csd,2))];

end


function csd = csdfromlfp_NicholsonFreeman(lfp,templfpData)

% Vaknin transform: (Vaknin 1985: add channels above and below assuming isoconductivity):
lfp_vaknin      = [lfp(1,:); lfp; lfp(end,:);];

% Filter temporally and spatially: is key to obtaining a smooth csd profile (filter over 3 channels and 5 time points)
spat_filter     = fspecial('gaussian',[3 5],1.1);
lfp_vaknin      = conv2(lfp_vaknin,spat_filter,'same');

% compute the CSD map:
csd             = CSD(lfp_vaknin',templfpData.fs(1),50*1e-6,'conductivity',0.4);
csd             = csd'; %Transpose

% Remove first and last 2 channels of csd map to ensure same size as original
% lfp matrix and set outer channels to zero (nonsensical csd values)
csd             = [zeros(1,size(csd,2));csd(3:end-2,:);zeros(1,size(csd,2))];

end

%%
%% Panel 1: raw traces:
% subplot(1,3,1);
% cla;
% 
% % If there is a task, then take a random lick, extract lfp surrounding it
% % to investigate movement artefacts:
% dirfiles = dir(fullfile(params.sessionRootDir,'ChangeDetection*'));
% 
% if length(dirfiles)>1
%     dirfiles = dirfiles(1);
% end
% 
% if ~isempty(dirfiles)
%     eventfilename = fullfile(params.sessionRootDir,dirfiles.name,'Events.nev');
%     %Settings for getting events
%     FieldSelection      = [1 0 1 0 1]; %Fields to be read: TimeStamps, TTLs, EventStrings
%     ExtractHeader       = 0;
%     ExtractMode         = 1; %Extract all data, from beginning to end
%     [TimeStamps, TTLs, EventStrings] =   Nlx2MatEV(eventfilename,FieldSelection,ExtractHeader,ExtractMode);
%     TTL_leftLick        = 1; %B0001xxxx
%     TTL_rightLick       = 2; %B0010xxxx
%     Port_leftLick      	= 2;
%     Port_rightLick      = 2;
%     EventStrings = EventStrings';
%     EventPort = NaN(1,length(EventStrings));
%     for i = 1:length(EventStrings)
%         EventPort(i) = str2double(EventStrings{i}(strfind(EventStrings{i},'port ')+5));
%     end
%     TS_Licks            = TimeStamps((EventPort == Port_leftLick & TTLs == TTL_leftLick) | (EventPort == Port_rightLick & TTLs == TTL_rightLick));
%     ChosenLick          = TS_Licks(ceil(length(TS_Licks)/randi(9)));
%     TimeWindow          = [ChosenLick + params.t_pre ChosenLick + params.t_post];
%     
%     %Parameters for filters and resampling:
%     prm.UseButter               = 0;
%     prm.UseKaiser               = 0;
%     prm.UseResampling           = 0;    %Resample data to lower frequency
%     %Extract the lfp in a subfunction:
%     templfpData_lick                = MOL_extractLFP(fullfile(params.sessionRootDir,'RawData'),params.OrigChannelSel,TimeWindow,prm);
%     
%     if params.flipsignLFP
%         templfpData_lick.signal = cellfun(@(x) -x, templfpData_lick.signal, 'UniformOutput',false);
%     end
%     
%     templfpData_lick.ts             = repmat(templfpData_lick.t_start(1):1/templfpData_lick.fs(1)*1e6:templfpData_lick.t_end(1),1,1);
%     
%     timeselec                   = templfpData_lick.ts> TimeWindow(1) & templfpData_lick.ts<TimeWindow(2);
%     offsetmat                   = repmat(templfpData.sortedChannelDepth(templfpData.sortedisgood),1,length(params.xtime));
% 
%     templfpData_lick.sortedSignal                = cell2mat(templfpData_lick.signal);
%     templfpData_lick.sortedSignal                = templfpData_lick.sortedSignal(templfpData.sortedChannelidx,:);
%     
%     %Plot:
%     baseline = repmat(mean(templfpData_lick.sortedSignal(templfpData.sortedisgood,timeselec),2),1,sum(timeselec));
%     handles = plot(params.xtime,(templfpData_lick.sortedSignal(templfpData.sortedisgood,timeselec)-baseline)*20e4 + offsetmat,'LineWidth',0.2); hold on;
%     %add legend to identify channels by their colors:
%     legend(flipud(handles),num2cell(num2str(flipud(templfpData.sortedChannelNum(templfpData.sortedisgood))),2),'Location','NorthEast','FontSize',350/params.nChannels);
% else
%     iTrial                      = randi(length(trialData.(params.AlignOnEvent)));
%     timeselec                   = templfpData.ts>trialData.(params.AlignOnEvent)(iTrial)+params.t_pre & templfpData.ts<trialData.(params.AlignOnEvent)(iTrial)+params.t_post;
%     offsetmat                   = repmat(templfpData.sortedChannelDepth(templfpData.sortedisgood),1,length(params.xtime));
%     
%     %Plot:
%     baseline = repmat(mean(templfpData.sortedSignal(templfpData.sortedisgood,timeselec),2),1,sum(timeselec));
%     handles = plot(params.xtime,(templfpData.sortedSignal(templfpData.sortedisgood,timeselec)-baseline)*20e4 + offsetmat,'LineWidth',0.2); hold on;
%     %add legend to identify channels by their colors:
%     legend(flipud(handles),num2cell(num2str(flipud(templfpData.sortedChannelNum(templfpData.sortedisgood))),2),'Location','NorthEast','FontSize',350/params.nChannels);
% end
% 
% %Make up:
% set(gca,'XTick',params.xticks,'XTickLabels',params.xticklabels)
% xlabel('Time from stimulus (ms)','FontSize', 15);
% ylabel('Channel','FontSize', 15)
% plot([0 0],params.Ylim_all,'k','LineWidth',2);
% xlim([params.xtime(1) params.xtime(end)]);
% title('Raw Traces','FontSize',15)
% ylabel('Channel Depth (in um from dura)','FontSize', 15)

% %%  MUA power:
% subplot(1,3,2);
% cla; hold all;
% 
% %Get power estimate for different frequency bands:
% %                         [pwr_out,pwr_f] = pwelch(sortedlfpmat,[],[],[],templfpData.fs(1,1),'power');
% [pwr_out,pwr_f]         = pwelch(templfpData.sortedSignal(:,1:20000)',[],[],[],templfpData.fs(1,1),'power');
% 
% hf_pwr                  = sum(pwr_out(pwr_f>params.hf_lower & pwr_f<params.hf_upper,:),1);
% hf_pwr                  = hf_pwr/max(hf_pwr(templfpData.sortedisgood)); %normalize to maximum
% plot(hf_pwr(templfpData.sortedisgood),templfpData.sortedChannelDepth(templfpData.sortedisgood),'-.k*','LineWidth',3)
% 
% linenoise_pwr                  = sum(pwr_out(pwr_f>49 & pwr_f<51,:),1);
% linenoise_pwr                  = linenoise_pwr/max(linenoise_pwr(templfpData.sortedisgood)); %normalize to maximum
% plot(linenoise_pwr(templfpData.sortedisgood),templfpData.sortedChannelDepth(templfpData.sortedisgood),'-.r','LineWidth',1)
% 
% theta_pwr                  = sum(pwr_out(pwr_f>8 & pwr_f<12,:),1);
% theta_pwr                  = theta_pwr/max(theta_pwr(templfpData.sortedisgood)); %normalize to maximum
% plot(theta_pwr(templfpData.sortedisgood),templfpData.sortedChannelDepth(templfpData.sortedisgood),'-.g','LineWidth',1)
% 
% xlabel('Normalized MUA power')
% legend({'MUA' '50 Hz' 'Theta'})
%%