function [] = boxplotPerc(data, percentile, varargin)

%==========================================================================
% FILE:     boxplotPerc.m
% AUTHOR:   Douglas Cook
% DATE:     3/12/21
% 
% PURPOSE:  A boxplot function which allows users to specify the percentile 
% values of boxplot whiskers while still using the standard MATLAB boxplot options. 
%
% INPUTS:   Use syntax exactly as with the standard "boxplot" command, but
% with a scalar or vector of percentile values as the second input.
%
%           data - a data matrix as with normal boxplot
%
%           percentile - a vector of 2 percentile values [lower, upper],
%           specified as percentages --OR-- a single scalar indicating the
%           desired span of percentiles. If the scalar is used, the program
%           solves for the necessary [lower, upper] values. 
%           Examples: [5 95] or equivalently [90]
%                     [2.5 97.5] or equivalently [95];
%
%           varargin - all the remaining boxplot options, according to the
%           typical MATLAB syntax. 
%
% OUTPUT:   a boxplot with the whiskers set to the specified percentile or
%           percentile range values.
% 
%
% EXAMPLE:  Here is a highly asymmetric data set along with a call to the boxplotPerc function:
%
%           data = normrnd(-1,1,100,5).^3;
%           boxplotPerc(data,90, 'notch','on','symbol','o','Widths',[0.1:0.1:0.5])
%           
%           Note: Visualization code showing the differences between standard
%           and percentage whiskers is shown below the main body.
%
% VERSION HISTORY
% V1.1 - DC: 3/12/21 - Adapted from a MATLAB file exchange file: https://www.mathworks.com/matlabcentral/fileexchange/22526-box-plot-with-whiskers-plotted-at-fixed-percentiles 
%==========================================================================


% ============= PART 1 - CREATE BOXPLOT ===================================
% This part deconstructs the boxplot options, forming them into a character arrar that can
% be re-evaluated using the same syntax as the original boxplot commmand
boxPlotCommand = ['boxplot(data'];          % initial part of the boxplot command
for i = 1:length(varargin)      
    if ischar(varargin{i})                  % if it's a character array, just insert with some syntax
        boxPlotCommand = [boxPlotCommand, ', ''', varargin{i}, ''''];
    elseif isa(varargin{i},'numeric')       % if it's a numeric array, brackets are needed
        boxPlotCommand = [boxPlotCommand, ', [', num2str(varargin{i}), ']'];
    end
end
boxPlotCommand = [boxPlotCommand, ');'];    % close off the command
eval(boxPlotCommand)                        % perform the boxplot command
%==========================================================================




% ============= PART 2 - MODIFY BOXPLOT LINES =============================

% Check the size of percentile variable and expand if needed
if length(percentile) == 1
    percentile = (100 - percentile)/2;              %
    percentile = [percentile, 100-percentile];      % symmetric percentile vector which will cover the requested span
end


% get the necessary handles.
b = get(gca,'children');                            % get children of the axes, which is a boxplot "group" of lines
c = get(b,'children');                              % the children of b are individual line arrays
% UW = findobj(c,'Tag','Upper Whisker');              % get the UpperWhisker line array
% UAJ = findobj(c,'Tag','Upper Adjacent Value');      % the Upper Crossbar line array
% OUT = findobj(c,'Tag','Outliers');                  % upper outliers
% LW = findobj(c,'Tag','Lower Whisker');              % Lower whisker line array
% LAJ = findobj(c,'Tag','Lower Adjacent Value');      % lower crossbar
% 
UW = findobj('Tag','Upper Whisker');              % get the UpperWhisker line array
UAJ = findobj('Tag','Upper Adjacent Value');      % the Upper Crossbar line array
OUT = findobj('Tag','Outliers');                  % upper outliers
LW = findobj('Tag','Lower Whisker');              % Lower whisker line array
LAJ = findobj('Tag','Lower Adjacent Value');      % lower crossbar

[~, m] = size(data);                                 % size of data array    
P = prctile(data, percentile);                      % percentiles
if m == 1, P = P(:);
lowerWhisker = P(1,:);                              % location of lower whiskers
upperWhisker = P(2,:);                              % location of upper whiskers

for i = 1:m

    UW_lims = get(UW(m-i+1),'ydata');                        % current upper whisker points
    set(UW(m-i+1),'ydata',[UW_lims(1) upperWhisker(i)]);     % set to same lower, new upper
    set(UAJ(m - i + 1),'ydata',[1 1]*upperWhisker(i));       % move the upper crossbar
    
    LW_lims = get(LW(m - i + 1),'ydata');                    % current lower whisker points
    set(LW(m-i+1),'ydata',[lowerWhisker(i), LW_lims(2)]);    % set to new lower, same upper       
    set(LAJ(m - i + 1),'ydata',[1 1]*lowerWhisker(i));       % move the lower crossbar
    
    index = data(:,i) > upperWhisker(i) | ...
            data(:,i) < lowerWhisker(i);   
            
    outliers = sort(data(index,i));                         % sort the outliers
    set(OUT(m-i+1),'Ydata',outliers);                       % reset the outliers
    xvals = OUT(m-i+1).XData;                               % get the xvalues for the outliers
    nx = length(xvals);                                     % number of xvalues
    nOut = length(outliers);                                % number of outliers
    index = randsample(nx,nOut,true);                       % resample the xvalues with replacement
    xvals = xvals(index);                                   % create a new vector of x values
    set(OUT(m-i+1),'Xdata',xvals);                          % reset the xvalues
end
% =========================================================================

end




%% VISUALIZATION OF DIFFERENCES BETWEEN STANDARD AND PERCENTILE WHISKERS
% data = (normrnd(1,1,1000,1).^3);
% ylimits = [-10, 40];
% 
% subplot(1,3,1)
% histogram(data, 100, 'orientation','horizontal')
% title('Histogram')
% ylim(ylimits)
% 
% subplot(1,3,2), boxplot(data, 'notch','on','symbol','.','Jitter',0.2)
% title('Standard Boxplot')
% ylim(ylimits)
% 
% subplot(1,3,3), boxplotPerc(data,[95], 'notch','on','symbol','.','Jitter',0.2)
% title('Boxplot with Whiskers at 95% coverage')
% ylim(ylimits)
