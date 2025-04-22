function y = show_sterr_p(x)

[~,p] = ttest(x,0,'dim',2);

disp([nanmean(x,2) nanstd(x,0,2)./sqrt(sum(~isnan(x),2)) p]);

y = [nanmean(x,2) nanstd(x,0,2)./sqrt(sum(~isnan(x),2)) p];