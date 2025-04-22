%plot_force_errorbars(X,Y,startbar,dist,args) plots Y with errorbars 
%against X. Observations are taken as columns of the 2D matrix Y. 
%startbar and startdist define where the errorbars start and how often
%they are plotted. NaNs are ignored.
%
%args contain plot properties and work as in plot.
%
%plot_force_errorbars(Y,startbar,dist,args) plots Y with errorbars against
%its index.
%
%Alkis
%5/23/2011

function plot_force_errorbars_dot_ci(X,Y,startbar,dist,varargin)

if ~isvector(X)
    varargin = {dist varargin{:}};
    dist = startbar;
    startbar = Y;
    Y=X;
    d = size(Y);
    X = 1:d(1);
end

X = X(:);

d = sum(~isnan(Y),2);
meanD = nanmean(Y,2);
errD = nanstd(Y,0,2)./sqrt(d);

X_ = X(startbar:dist:end);
meanD_ = meanD(startbar:dist:end);
errD_ = errD(startbar:dist:end);

hold on;
plot(X,meanD,'.','markersize',20,varargin{:});
plot([X_ X_]',[meanD_+errD_ meanD_-errD_]',varargin{:});