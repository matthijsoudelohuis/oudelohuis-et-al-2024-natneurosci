function MOL_identifyMovArtefac(sessionData)

global mainfig params temptrialData templfpData

nSessions = length(sessionData.session_ID);

for iSes = 1:nSessions
    fprintf('Computing CSD for session %d/%d\n',iSes,nSessions)
    
    if ~(params.overwrite==0 && exist(fullfile(params.savedir,['MovArt_' sessionData.session_ID{iSes} '.mat']),'file'))
        
        %Get the relevant data for each session individually with the actual lfp data:
        [Data]              = MOL_GetData('E:','CHDET',{'ChangeDetectionConflict' 'ChangeDetectionConflictDecor' 'VisOnlyTwolevels'},{},sessionData.Rec_datetime(iSes),{'sessionData' 'trialData_newtrials' 'lfpData'});
        tempsessionData     = Data.sessionData;
        temptrialData       = Data.trialData;
        templfpData         = Data.lfpData;
        temptrialData       = MOL_RemoveLastnTrials(temptrialData,5);
        
        keepidx = strcmp(templfpData.area,params.area);
        lfpfields = fieldnames(templfpData);
        for iF = 1:length(lfpfields)
            templfpData.(lfpfields{iF}) = templfpData.(lfpfields{iF})(keepidx,:);
        end
        
        templfpData.sortedChannelDepth = -templfpData.sortedChannelDepth;
        if all(templfpData.sortedChannelDepth>0)
            templfpData.sortedChannelDepth = -templfpData.sortedChannelDepth;
        end
        
        if ~isempty(templfpData.signal) && any(templfpData.lfpgood)
            
            if params.flipsignLFP
                templfpData.signal = cellfun(@(x) -x, templfpData.signal, 'UniformOutput',false);
            end
            
            switch tempsessionData.(sprintf('Probe%d_Config',unique(templfpData.probeIdx))){1}
                case 'A4x8-5mm-100-200-177-CM32'
                    params.nChannels           = 32;
                    params.nShanks             = 4;
                    params.nChannelsShank      = 8;
                    params.intersitedistance   = 100;
                case 'A1x32-Poly2-5mm-50s-177-CM32'
                    params.nChannels           = 32;
                    params.nShanks             = 2;
                    params.nChannelsShank      = 16;
                    params.intersitedistance   = 50;
                case 'A1x64-Poly2-6mm-23s-160'
                    params.nChannels           = 64;
                    params.nShanks             = 2;
                    params.nChannelsShank      = 32;
                    params.intersitedistance   = 46;
            end
            
            params.Ylim_all                    	= [min(templfpData.sortedChannelDepth)-25 max(templfpData.sortedChannelDepth)+25];
            params.nChannels                    = length(templfpData.ch);
            
            templfpData.ts                      = repmat(templfpData.t_start(1):1/templfpData.fs(1)*1e6:templfpData.t_end(1),1,1);
            
            templfpData.sortedChannelNum = mod(templfpData.sortedChannelNum-1,params.nChannels)+1;
            templfpData.sortedChannelidx = mod(templfpData.sortedChannelidx-1,params.nChannels)+1;
            
            templfpData.sortedSignal            = cell2mat(templfpData.signal);
            templfpData.sortedSignal            = templfpData.sortedSignal(templfpData.sortedChannelidx,:);
            
            params.ChannelSel                   = 1:params.nChannels;
            
            mainfig     = figure; set(gcf,'units','normalized','Position',[0.1 0.2 0.84 0.62],'color','w');
            set(gcf,'defaultAxesFontSize',15)
            title(sprintf('CSD overview - %s - %s',tempsessionData.mousename{1},strrep(tempsessionData.Rec_datetime{1},'_',' ')));
            
            params.ChannelSel                  = 1:params.nChannels;
            
            templfpData.sortedisgood                = true(size(templfpData.ch));
            
            %% Panel 1, 2: CSD
            MOL_UpdateCSDfig();
            
            %% Movement artefact selection panel:
            
            params.goodbadcolors      = {[.9 .2 .3] [.3 .9 .3]};
            xmin                = 0.91;
            xmax                = 0.99;
            params.ymax         = 0.96;
            params.ymin         = 0.1;
            
            y_width             = (params.ymax-params.ymin)/params.nChannelsShank;
            x_width             = (xmax-xmin)/params.nShanks;
            
            for iCh = 1:params.nChannels
                x_pos           = (templfpData.ChannelX(iCh)-min(templfpData.ChannelX))/ (max(templfpData.ChannelX)-min(templfpData.ChannelX));
                x_pos           = xmin + x_pos*(params.nShanks-1)*x_width;
                y_pos           = (templfpData.ChannelY(iCh)-min(templfpData.ChannelY))/ (max(templfpData.ChannelY)-min(templfpData.ChannelY));
                y_pos           = params.ymin + (1-y_pos) * (params.ymax-params.ymin-y_width);
                
                uicontrol('Parent',mainfig,'Style','togglebutton','units','normalized',...
                    'Position',[x_pos y_pos x_width y_width],'String',sprintf('Ch %3.0f',iCh),...
                    'backgroundcolor',params.goodbadcolors{templfpData.sortedisgood(iCh)+1},...
                    'Value',1,'UserData',iCh,'Callback',@togglebutton_Callback);
            end
            
            
            %% Finalize button:
            x_pos = 0.1;
            callbackstr = strcat('global mainfig; mainfig.UserData = 0; uiresume(mainfig);');
            uicontrol('Parent',mainfig,'Style','push','units','normalized',...
                'Position',[0.84 0.01 0.08 0.06],'String','Next','fontsize',16,...
                'backgroundcolor',[0.2,0.6,0.9],'Callback',callbackstr); x_pos = x_pos+x_width*1.1;
            callbackstr = strcat('global mainfig; mainfig.UserData = 1; uiresume(mainfig);');
            uicontrol('Parent',mainfig,'Style','push','units','normalized',...
                'Position',[0.92 0.01 0.08 0.06],'String','All Bad','fontsize',16,...
                'backgroundcolor',[0.9,0.3,0.2],'Callback',callbackstr);
            
            %Wait for data button to be pressed: callback is continue with script
            uiwait(mainfig)
            
            %% Saving the data:
            if mainfig.isvalid
                if mainfig.UserData
                    templfpData.sortedisgood = zeros(size(templfpData.session_ID));
                end
                templfpData     = rmfield(templfpData,{'signal','sortedSignal'});
                
                list            = templfpData.channel_ID(templfpData.sortedChannelidx(templfpData.sortedisgood==1));
                
                save(fullfile(params.savedir,['MovArt_' sessionData.session_ID{iSes} '.mat']),'list','templfpData')
            end
            if mainfig.isvalid
                close(mainfig);
            end
            
        end
    end
end

end

function togglebutton_Callback(hObject,eventdata,handles) %#ok<INUSD>
% hObject    handle to togglebutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

button_state = get(hObject,'Value');
global templfpData params
templfpData.sortedisgood(templfpData.sortedChannelNum==hObject.UserData) = button_state;
set(hObject,'BackgroundColor',params.goodbadcolors{button_state+1});
MOL_UpdateCSDfig();
fprintf('Figures updated for channel %d\n',hObject.UserData)
end

