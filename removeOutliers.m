function [x, N, List] = removeOutliers(x,nIQR)
xIQR = (x - median(x,2))./ iqr(x,2);
N = sum(sum(xIQR>nIQR));
N = [N, N/numel(xIQR)];
List = [];
x( abs(xIQR)>nIQR ) = NaN;
end