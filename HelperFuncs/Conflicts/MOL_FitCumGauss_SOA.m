
function [mu, stddev, guessrate, lapserate, curve] = MOL_FitCumGauss_SOA(xdata,ydata)

% Function: F=@(g,l,u,v,x) g+(1-g-l)*0.5*(1+erf((x-u)/sqrt(2*v^2)));

% The psychometric function presented in Wichmann and Hill's paper
% is a cumulative gaussian function with 4 parameters:
% Mean (u): The mean value of the distribution representing subject bias.
% Standard deviation (v): The variation of the distribution representing
%   the subjects discrimination sensitivity.
% Guess rate (g) and lapse rate (l): Two additional parameters representing
%   the subjects fallibility (ie. potential inability to ever reach 100%
%   performance) at each end of the distribution/stimulus spectrum.

% Set limits and starting positions for fit:
UL =       [1,      1,    0.3,    0.2];     % Upper limits for g, l, u ,v
SP =       [0.2,    0.2,    0,      0.1];      % Start points for g, l, u ,v
LL =       [0,      0,      -0.3,   0];        % Lower limits for  g, l, u ,v

%Check range of data:
% if min(ydata)<0 || max(ydata)>1
%     % Attempt to normalise data to range 0 to 1
%     ydata = ydata/(mean(ydata)*2);
% end

% Prepare fitting function
% F=@(g,l,u,v,x) g+(1-g-l)*0.5*(1+erf((-x-u)/sqrt(2*v^2)));
% F=@(g,l,u,v,x) g+(1-g-l)*(1+erf((-x-u)/sqrt(2*v^2)));
% F=@(g,l,u,v,x) 1-g+(1-g-l)*-(1+erf((x-u)/sqrt(2*v^2)));
F=@(g,l,u,v,x) 1-g+(2-g-l)*0.5*-(1+erf((x-u)/sqrt(2*v^2)));

% SPs and limits specified, use while fitting
ffit=fit(xdata,ydata',F,'StartPoint',SP,'Upper',UL,'Lower',LL);

mu          = ffit.u;
stddev      = ffit.v;
guessrate   = ffit.g;
lapserate   = ffit.l;

% Create a new xAxis with higher resolution
% fineX = linspace(min(xdata),max(xdata),numel(xdata)*50);
fineX = linspace(min(xdata),max(xdata),1000);

% Generate curve from fit
curve = feval(ffit, fineX);
curve = [fineX', curve];

end



