function [csdmat,erpmat] = MOL_CSD_trials(params, events_ts, lfpData)

% function [csdmat,erpmat] = calcCSD(params,lfpData,trialData)
%make butter filter:
[params.B_butt,params.A_butt]   = butter(params.ord_butter,params.lp_butter/(lfpData.fs(1)/2));

if params.UseKaiser
    [N,Wn,beta,ftype]   = kaiserord([params.hp_kaiser-params.dp_kaiser params.hp_kaiser params.lp_kaiser params.lp_kaiser+params.dp_kaiser],[0 1 0],[1e-6 1 1e-6],lfpData.fs(1));
    B_kai               = fir1(N,Wn,ftype,kaiser(N+1,beta),'noscale');
    A_kai               = 1;
end

nEvents                         = length(events_ts);
erpmat                          = zeros(params.nChannels,length(params.xtime),nEvents); %init lfp matrix

%Get lfp per trial and apply low-pass filter:
for ev = 1:nEvents
    
    tempsignal      = lfpData.sortedSignal(:,lfpData.ts>events_ts(ev)+params.t_pre & lfpData.ts<events_ts(ev)+params.t_post);
    erpmat(:,:,ev)  = tempsignal(:,1:length(params.xtime));
    
    if params.UseButter
        for iCh = 1:params.nChannels
            erpmat(iCh,:,ev)            = filtfilt(params.B_butt,params.A_butt,erpmat(iCh,:,ev));
        end
    end
    
    if params.UseKaiser
        for iCh = 1:params.nChannels
            erpmat(iCh,:,ev)            = filtfilt(B_kai,A_kai,erpmat(iCh,:,ev));
        end
    end
end

% Replace bad channels with mean of adjacent channels
for iCh = find(~lfpData.sortedisgood)'
    chbelow          = find(lfpData.sortedChannelDepth > lfpData.sortedChannelDepth(iCh) & lfpData.sortedisgood,1,'first');
    chabove          = find(lfpData.sortedChannelDepth < lfpData.sortedChannelDepth(iCh) & lfpData.sortedisgood,1,'last');
    erpmat(iCh,:,:)  = mean(erpmat([chbelow chabove],:,:),1);
end

%Loop over the different shanks (CSD does not work well with
%considering no horizontal spacing)
for iShank = 1:params.nShanks
    ShankChannelSel                     = params.ChannelSel([1:params.nChannelsShank] + params.nChannelsShank*(iShank-1));
    
    idx                                 = ismember(lfpData.sortedChannelNum,ShankChannelSel);
%     sortedShankChannelNum               = lfpData.sortedChannelNum(idx);
    sortedShankChannelDepth             = lfpData.sortedChannelDepth(idx);
    %     sortedShankChannelIdx       = lfpData.sortedChannelNum(ismember(lfpData.sortedChannelNum,ShankChannelSel));
    
    shanklfpmat                         = erpmat(idx,:,:);

%     switch params.csdmethod
%         case 'DoubleDiff'
%             csdmat                             = csdfromlfp_doublediff(mean(shanklfpmat,3));
%         case 'NicholsonFreeman'
%             csdmat                             = csdfromlfp_NicholsonFreeman(mean(shanklfpmat,3),lfpData);
%     end
    
    for ev = 1:nEvents
        switch params.csdmethod
            case 'DoubleDiff'
%                 csdmat                             = csdfromlfp_doublediff(mean(shanklfpmat,3));
            case 'NicholsonFreeman'
                csdmat(:,:,ev)                     = csdfromlfp_NicholsonFreeman(shanklfpmat(:,:,ev),lfpData);
        end
    end
    ChannelY_shank{iShank}             = sortedShankChannelDepth;
    csdmat_shank{iShank}               = csdmat;
end

%Baseline subtraction of csd / erp:
for iShank = 1:params.nShanks
    csdmat_shank{iShank}               = csdmat_shank{iShank} - repmat(mean(mean(csdmat_shank{iShank}(:,params.xtime<-0.025),2),3),1,length(params.xtime),nEvents);
end

% baseline              = repmat(mean(erpmat(:,params.xtime<0),2),1,length(params.xtime),1);
% erpmat                = erpmat - baseline;
baseline               = repmat(mean(mean(erpmat(:,params.xtime<0),2),3),1,length(params.xtime),nEvents);
erpmat                = erpmat - baseline;

%interpolate two shanks to same channels and then average across shanks:
for ev = 1:nEvents
    csdmat_shank_int{1}(:,:,ev)        = tointerpol2(csdmat_shank{1}(:,:,ev),ChannelY_shank{1}',lfpData.sortedChannelDepth);
    csdmat_shank_int{2}(:,:,ev)        = tointerpol2(csdmat_shank{2}(:,:,ev),ChannelY_shank{2}',lfpData.sortedChannelDepth);
end
csdmat                 = cat(4,csdmat_shank_int{1},csdmat_shank_int{2});
csdmat                 = nanmean(csdmat,4);  

% csdmat_shank{1}        = csdmat_shank{1}/max(max(csdmat_shank{1}));
% csdmat_shank{2}        = csdmat_shank{2}/max(max(csdmat_shank{2}));

% csdmat                 = cat(3,csdmat_shank{1},csdmat_shank{2});
% csdmat                 = nanmean(csdmat,3);

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


function csd = csdfromlfp_NicholsonFreeman(lfp,lfpData)

% Vaknin transform: (Vaknin 1985: add channels above and below assuming isoconductivity):
lfp_vaknin      = [lfp(1,:); lfp; lfp(end,:);];

% Filter temporally and spatially: is key to obtaining a smooth csd profile (filter over 3 channels and 5 time points)
spat_filter     = fspecial('gaussian',[3 5],1.1);
lfp_vaknin      = conv2(lfp_vaknin,spat_filter,'same');

% compute the CSD map:
csd             = CSD(lfp_vaknin',lfpData.fs(1),50*1e-6,'conductivity',0.4);
csd             = csd'; %Transpose

% Remove first and last 2 channels of csd map to ensure same size as original
% lfp matrix and set outer channels to zero (nonsensical csd values)
csd             = [zeros(1,size(csd,2));csd(3:end-2,:);zeros(1,size(csd,2))];

end




% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% %% Get input arguments:
% sessionData     = varargin{1};
% trialData       = varargin{2};
% % spikeData       = varargin{3};
% lfpData         = varargin{4};
% 
% 
% %% DFT filter on 50 Hz:
% % L = size(lfpData.signal,2);
% % f = lfpData.fs*(0:(L/2))/L;
% % for iCh = 1:length(lfpData.ts)
% %     fprintf('Filtering 50Hz noise channel %d\n',iCh)
% %     h = fft(lfpData.signal(iCh,:));
% %     h(find(f>49 & f<51)) = complex(zeros(numel(find(f>49 & f<51)),1));
% %     lfpData.signal(iCh,:) = ifft(h);
% % end
% 
% %%
% nPreviousChannels = 0; %counter that keeps track of cumulative channel number
% for iProbe = 1:3
%     if sessionData.(sprintf('Probe%d',iProbe))
%         switch sessionData.(sprintf('Probe%d_Config',iProbe))
%             case 'A1x32-Poly2-5mm-50s-177-CM32'
%                 nChannels           = 32;
%                 nShanks             = 2;
%                 nChannelsShank      = 16;
%                 intersitedistance   = 50;
%             case 'A4x8-5mm-100-200-177-CM32'
%                 nChannels           = 32;
%                 nShanks             = 4;
%                 nChannelsShank      = 8;
%                 intersitedistance   = 100;
%             case 'A1x64-Poly2-6mm-23s-160'
%                 nChannels           = 64;
%                 nShanks             = 2;
%                 nChannelsShank      = 32;
%                 intersitedistance   = 46;
%         end
%         
%         ChannelSel                      = [1:nChannels] + nPreviousChannels;
%         nPreviousChannels               = nPreviousChannels + nChannels;
%         xtime                           = (params.t_pre:1e6/lfpData.fs:params.t_post-0.01e6) * 1e-6;
%         nEvents                         = length(events_ts);
%         lfpmat                          = zeros(nChannels,length(xtime),nEvents); %init lfp matrix
%         
%         %%Get lfp per trial
%         %make butter filter:
%         [params.B_butt,params.A_butt]     = butter(params.ord_butter,params.lp_butter/(lfpData.fs(1)/2));
% 
%         for ev = 1:nEvents
%             tempsignal      = lfpData.signal(ChannelSel,lfpData.ts{1,:}>events_ts(ev)+params.t_pre & lfpData.ts{1,:}<events_ts(ev)+params.t_post);
%             lfpmat(:,:,ev)  = tempsignal(:,1:length(xtime));
%             
%             if params.UseButter
%                 for iCh = 1:nChannels
%                     lfpmat(iCh,:,ev)            = filtfilt(params.B_butt,params.A_butt,lfpmat(iCh,:,ev));
%                 end
%             end
%         end
%         
%         %Loop over the different shanks (CSD does not work well with
%         %considering no horizontal spacing)
%         for iShank = 1:nShanks
%             ShankChannelSel                     = ChannelSel([1:nChannelsShank] + nChannelsShank*(iShank-1));
%             
%             [ShankChannelY,sortidx]             = sort(-lfpData.ChannelY(ShankChannelSel),'ascend');
%             sortidx                             = sortidx + nChannelsShank*(iShank-1);
%             
%             shanklfpmat                         = NaN(nChannelsShank,size(lfpmat,2),size(lfpmat,3));
%             shanklfpmat(:,:,:)                  = lfpmat(sortidx,:,:);
%             
%             %Init csd matrix:
%             csd                                 = zeros(size(shanklfpmat));
%             for ev = 1:size(shanklfpmat,3)
%                 csd(:,:,ev)                     = csdfromlfp(shanklfpmat(:,:,ev));
%             end
%             csdmat                             = mean(csd,3);
%             csdmat                             = csdmat / (intersitedistance/100);
% 
%             %Baseline subtraction of csd:
%             csdmat_shank{iShank}               = csdmat - repmat(mean(csdmat(:,xtime<-0.025),2),1,length(xtime));
%             ChannelY_shank{iShank}              = ShankChannelY;
%         end
%         
%         %Interpolate values for visualization purposes:
%         newChannelY = linspace(min(min([ChannelY_shank{:}])),max(max([ChannelY_shank{:}])),100);
% 
%         if params.interpolate
%             csdmat_shank{1}                            = tointerpol2(csdmat_shank{1},ChannelY_shank{1}',newChannelY);
%             csdmat_shank{2}                            = tointerpol2(csdmat_shank{2},ChannelY_shank{2}',newChannelY);
%         end
%         
%         csdmat = cat(3,csdmat_shank{1},csdmat_shank{2});
%         csdmat = nanmean(csdmat,3);
%         
%         params.cscale = [-max(max(csdmat))*0.9 max(max(csdmat))*0.9];
%         
%         [ChannelY_meanLFP,sortidx]  = sort(-lfpData.ChannelY(ChannelSel));
%         lfpmat                      = lfpmat(sortidx,:,:);
%         
%         csdfig = figure; set(gcf,'units','normalized','Position',[0.05 0.4 0.9 0.4],'color','w');
%         set(gcf,'defaultAxesFontSize',15)
%         suptitle(sprintf('CSD %s - %s',sessionData.mousename,sessionData.(sprintf('Probe%d_Area',iProbe))));
%         
%         subplot(1,3,1)
%         offsetmat = repmat(ChannelY_meanLFP,1,size(lfpmat(:,:,1),2));
%         plot(xtime,mean(lfpmat(:,:,:),3)*20e4 + offsetmat); hold on;
%         xlabel('Time from stimulus (s)','FontSize', 15)
%         ylabel('Channel','FontSize', 15)
%         plot([0 0],ylim,'k','LineWidth',2);
%         xlim([xtime(1) xtime(end)]);
% %         ylim([min(ChannelY_meanLFP)-50 max(ChannelY_meanLFP)+50]);
%         title('Mean LFP (<50Hz) to Checker reversal','FontSize',15)
%         ylabel('Channel Depth (in um from dura)','FontSize', 15)
%         
%         subplot(1,3,2)
%         set(gcf,'defaultAxesFontSize',15)
% %         imagesc(xtime,flipud(newChannelY),flipud(csdmat),params.cscale)
%         imagesc(xtime,newChannelY,csdmat,params.cscale)
%         xlabel('Time from stimulus (s)','FontSize', 15)
%         set(gca,'YDir','normal')
% 
%         switch params.colormap
%             case 'redblue'
%                 h       = coolwarm + repmat(0.1-abs(linspace(-0.1,.1,64))',1,3);
%                 h(h>1)  = 1;
%                 if params.sinkhot
%                     h = flipud(h);
%                 end
%                 colormap(h);
%             case 'parula'
%                 if params.sinkhot
%                     colormap(flipud(parula));
%                 else                     colormap(parula);
%                 end
%         end
%         ticks = linspace(params.cscale(1),params.cscale(2),9);
%         c = colorbar('Ticks',ticks,'TickLabels',num2cell(ticks*1e6),'FontSize',10);
%         c.Label.String = 'Sink                    (uA/mm3)                 Source';
%         
%         %Get power estimate for different frequency bands:
%         [pwr_out,pwr_f] = pwelch(lfpData.signal(ChannelSel,:)',[],[],[],lfpData.fs(1,1),'power');
%         
%         hf_pwr = sum(pwr_out(pwr_f>500 & pwr_f<5000,:),1);
%         hf_pwr = hf_pwr/max(hf_pwr); %normalize to maximum
%         
%         subplot(1,3,3)
%         [ChannelY_meanLFP,sortidx]  = sort(-lfpData.ChannelY(ChannelSel));
%         hf_pwr = hf_pwr(sortidx);
%         plot(hf_pwr,ChannelY_meanLFP)
% 
%         for i = 1:3
%             subplot(1,3,i)
%             ylim([min(ChannelY_meanLFP)-50 max(ChannelY_meanLFP)+50]);
%             
%             Vert_Ticks = -2000:100:400;
%             Vert_Ticks = Vert_Ticks(Vert_Ticks>min(ChannelY_meanLFP)-50 & Vert_Ticks<max(ChannelY_meanLFP)+50);
%             
%             set(gca, 'YTick', Vert_Ticks, 'YTickLabels', Vert_Ticks,'FontSize', 15)
%         end
%                 
%     end
% end
% 
% end
% 
% 
% function csd = csdfromlfp(lfp)
% 
% % Vaknin transform: (Vaknin 1985: add channels above and below assuming isoconductivity):
% lfp_vaknin      = [lfp(1,:); lfp(1,:); lfp; lfp(end,:); lfp(end,:)];
% 
% % Filter temporally and spatially: is key to this
% spat_filter     = fspecial('gaussian',[3 5],1.1);
% lfp_vaknin      = conv2(lfp_vaknin,spat_filter,'same');
% 
% % compute the CSD map
% csd             = diff(lfp_vaknin,2,1); %Taking the double derivative from up to down
% 
% % Remove first and last channels of csd map to ensure same size as original
% % lfp matrix and set outer channels to zero (nonsensical csd values)
% csd             = [zeros(1,size(csd,2));csd(3:end-2,:);zeros(1,size(csd,2))];
% 
% end