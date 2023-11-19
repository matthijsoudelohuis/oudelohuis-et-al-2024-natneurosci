function [data_highpercent,data_highrank,data_lowpercent,data_lowrank] = boxplot_percentwhisk(data,percen,criterion,group,group_list,column,fignum,xlabels,ylabels,titles)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% boxplot_percenwhisk.m
%
% Robert Elleman 12/24/2008
%
%
% This subroutine plots a boxplot figure but plots the whiskers according
% to a specified percentile not according to the way matlab does it by
% default. The outliers are only displayed outside the whisker.
%
% Inputs:
%        1. data:  
%               Input data. This can be a matrix. The first column is
%               the grouping dimension (like time).
%        2. percen:
%               Percentile for plotting the whisker. You should input the
%               upper percentile as 0-100. For example, input "98" if you
%               want the lower whisker to go to the 2nd percentile and the
%               upper whisker to go to the 98th percentile
%        3. criterion:
%               Value you would like to plot as a dotted line. For example,
%               if it is bad for data to be above 150, then "150" as
%               criterion will plot a dotted line at 150 from left to right
%               on your plot. If you don't want this, you should comment
%               out that part below.
%        4. group:
%               A Matlab group that segregates the input data into
%               groupings for separate boxes on the boxplot.
%        5. group_list:
%               The list of groups for different boxes in the boxplot 
%        6. column:
%               Column in data to plot
%        7. fignum:
%               Figure number for plotting
%        8. xlabels:
%               Text for xlabels in figure
%        9. ylables:
%               Text for ylabels in figure
%        10. titles:
%               Text for title
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% first plot out the Matlab version of boxplot so that you can get the
% handles.
figure(fignum)
boxplot(data(:,column),group,'medianstyle','target','notch','on','symbol','rx','outliersize',8)


% get the necessary handles.
b = get(gca,'children'); 
c = get(b,'children'); 
UW = findobj(c,'Tag','Upper Whisker');
UAJ = findobj(c,'Tag','Upper Adjacent Value');
OUT = findobj(c,'Tag','Outliers');
LW = findobj(c,'Tag','Lower Whisker');
LAJ = findobj(c,'Tag','Lower Adjacent Value');

% for each grouping along the x-axis, replot the whisker to be based on
% percentiles of the data. Also, only the outliers outside the new whiskers
% are plotted. For the percentile, it is calculated and then rounded up to
% the next nearest data value. For example, if the 98th percentile is
% between the 4th highest and 5th highest data value, then the percentile
% is set to be the 4th highest data value.

%keyboard

for i=1:length(group_list)
    group_index = ismember(group,group_list(i));
    data_sub = data(group_index,column);    
    [data_highpercent(i) data_highrank(i)] = percentile(data_sub,percen);
    [data_lowpercent(i) data_lowrank(i)] = percentile(data_sub,100-percen);
    UW_lims = get(UW(length(group_list)-i+1),'ydata');
    set(UW(length(group_list)-i+1),'ydata',[UW_lims(1) data_highpercent(i)])
    set(UAJ(length(group_list)-i+1),'ydata',[data_highpercent(i) data_highpercent(i)])
    LW_lims = get(LW(length(group_list)-i+1),'ydata');
    set(LW(length(group_list)-i+1),'ydata',[data_lowpercent(i) LW_lims(2)])
    set(LAJ(length(group_list)-i+1),'ydata',[data_lowpercent(i) data_lowpercent(i)])
    outliers = get(OUT(length(group_list)-i+1),'ydata');
    outlier_index = find(outliers == data_highpercent(i));
    outliers(1:outlier_index) = data_highpercent(i);
    set(OUT(length(group_list)-i+1),'ydata',outliers);
end

% Here a horizontal dotted line is plotted to mark the criterion. Comment
% if you are not interested in this feature.
hold on
plot(0:length(group_list)+1,repmat(criterion,[length(group_list)+2 1]),'k--','MarkerSize',8)
hold off 

% make labels and the title
xlabel(xlabels,'FontSize',12,'FontWeight','bold')
ylabel(ylabels,'FontSize',12,'FontWeight','bold')
title(titles,'FontSize',14,'FontWeight','bold')


function [out,num] = percentile(data,kpercent)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Percentile Calculation Example
% inputs:
%        1. data: input dataset expected as a column
%        2. kpercent: single value for percentage [0-100]
% outputs:
%        1. percentile "kpercent" from data 
%        2. rank of output value. Like if it is the 4th highest, then
%           this value is "4"
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% STEP 1 - rank the data
y = sort(nanfree(data,1));

% STEP 2 - find k% (k /100) of the sample size, n.
k = kpercent/100;
result = k*length(y);

% STEP 3 - take the value higher than kpercent

if (~isnan(y))
    out = y(ceil(result));
    num = length(y) - floor(result);
else
    out = NaN;
    num = NaN;
end


function out = nanfree(in,column)

%
% This function strips out the rows that have NaN in the specified
% column(s)
%

out = in;
if isempty(column)
    echo "please specify columns"
elseif (length(column) >= 1)
    index = isnan(mean(in(:,column),2));
    out(index,:) = [];
else
     echo "please specify positive columns"   
end



