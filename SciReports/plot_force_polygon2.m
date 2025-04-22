function plot_force_polygon2(x,y,color)

%d = size(y);
d = sum(~isnan(y),2); %to make this more accurate
ymean = nanmean(y,2);
yerr = nanstd(y,0,2)./sqrt(d);

interruptions = [1; find(isnan(ymean)); length(ymean)+1];

if isempty(find(isnan(ymean)))
    fill([x flip(x)]',[ymean-yerr; flip(ymean+yerr)],color,'EdgeColor',color);
else
    for i = 1:length(interruptions)-1
        t_i =  (interruptions(i)+1):(interruptions(i+1)-1);
        fill([x(t_i) flip(x(t_i))]',[ymean(t_i)-yerr(t_i); flip(ymean(t_i)+yerr(t_i))],color,'EdgeColor',color);
    end
end