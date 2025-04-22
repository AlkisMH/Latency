%plot_errorbars(X,Y,bartype,bartypearg,varargin) plots Y with errorbars
%against X. Observations are taken as columns of the 2D matrix Y.
%Add more explanation
% NaNs are ignored.
%
%args contain plot properties and work as in plot.
%
%plot_force_errorbars(Y,bartype,bartypearg,varargin) plots Y with errorbars against
%its index.
%
%Alkis
%Updated 4/7/2025

function plot_errorbars(X,Y,bartype,bartypearg,varargin)

if ~isvector(X)
    varargin = {bartypearg varargin{:}};
    bartypearg = bartype;
    bartype = Y;
    Y=X;
    d = size(Y);
    X = 1:d(1);
end

X = X(:);

d = sum(~isnan(Y),2);
meanD = nanmean(Y,2);
switch bartype
    case 'sem'
        errD = nanstd(Y,0,2)./sqrt(d);
    case 'std'
        errD = nanstd(Y,0,2);
    case 'ci'
        %Confidence interval based on SEM
        if isempty(bartypearg)
            bartypearg = 95;
        end
        errD = icdf('normal',1-(1-bartypearg/100)/2,0,1)*nanstd(Y,0,2)./sqrt(d);
end
hold on;
plot(X,meanD,'.','markersize',20,varargin{:});
plot([X X]',[meanD+errD meanD-errD]',varargin{:});