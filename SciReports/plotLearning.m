function x = plotLearning(z, nIQR, c, window)

dList = [0 85 300];
dList2 = [25 85 300];

%Change, 3/5/2024: do outlier rejection on total and explicit instead of
%implicit and explicit, and then calculate implicit

xo  = z.x.dat.tot .* sign(z.x.vmr);
xo0 = z.x.dat.tot_base .* sign(z.x.vmr);

%Outlier removal - do separately for each condition
xe  = z.x.dat.exp .* sign(z.x.vmr);
xe0 = z.x.dat.exp_base .* sign(z.x.vmr);

for m = 1:3
    i_m = find(z.x.delay==dList(m));
    xe(:,i_m) = removeOutliers(xe(:,i_m),nIQR);
    xe0(:,i_m) = removeOutliers(xe0(:,i_m),nIQR);
    xo(:,i_m) = removeOutliers(xo(:,i_m),nIQR);
    xo0(:,i_m) = removeOutliers(xo0(:,i_m),nIQR);
end
%x.NNi = sum(isnan([xi0;xi]));
x.NNe = sum(isnan([xe0;xe]));
x.NNo = sum(isnan([xo0;xo]));

%Both combined
x.NNr = sum(isnan([xe0;xe].*[xo0;xo]));

%Baseline substraction
xee = [xe0; xe] - nanmean(xe0(2:16,:));
xoo = [xo0; xo] - nanmean(xo0(2:16,:));

xii = xoo-xee;

i_break = [1 27 57 87 117 137 167];

%To compare (non-flipped) baselines
x.o{1}.base_signed = nanmean(xo0(2:16,:)) .* sign(z.x.vmr);
x.i{1}.base_signed = (nanmean(xo0(2:16,:))-nanmean(xe0(2:16,:))) .* sign(z.x.vmr);
x.e{1}.base_signed = nanmean(xe0(2:16,:)) .* sign(z.x.vmr);

%% Main adaptation plot

tt = (1:196)-16;
xii_nb = xii; xii_nb(i_break,:) = NaN;
xee_nb = xee; xee_nb(i_break,:) = NaN;

lc.i = xii_nb;
lc.e = xee_nb;
lc.o = xii_nb+xee_nb;

conds_ = 'oie';
cols_ = {'purps','blues','reds'};

figure(101);
set(101,'Position',[200,200,800,200]);

for n = 1:3
    subplot(1,3,n);hold on;
    plot(tt,tt*0,'k--');
    for k = 1:length(i_break)
        plot(tt([i_break(k) i_break(k)]),[-5 36],'k--')
    end
    for m = 1:3
        plot_force_polygon2(tt,lc.(conds_(n))(:,z.x.delay==dList(m)),0.5*c.(cols_{n})(m,:)+0.5);
    end
    for m = 1:3 %This is a separate loop so all averages are on top
        plot(tt,nanmean(lc.(conds_(n))(:,z.x.delay==dList(m)),2),'color',c.(cols_{n})(m,:),'linewidth',1);
    end
    xlim(tt([1 end]));ylim([-4 36]);xlim([tt(1) 120]);
    set(gca,'Xtick',0:20:120);
end

%Relearning figure
i_break_r = [1 31];
tt_r = 1:60;

figure(102);
set(102,'Position',[200,200,800,200]);

for n = 1:3
    subplot(1,3,n);hold on;
    for k = 1:length(i_break_r)
        plot(tt_r([i_break_r(k) i_break_r(k)]),[-5 36],'k--')
    end
    for m = 1:3
        plot_force_polygon2(tt_r,lc.(conds_(n))(137:end,z.x.delay==dList(m)),0.5*c.(cols_{n})(m,:)+0.5);
    end
    for m = 1:3 %This is a separate loop so all averages are on top
        plot(tt_r,nanmean(lc.(conds_(n))(137:end,z.x.delay==dList(m)),2),'color',c.(cols_{n})(m,:),'linewidth',1);
    end
    xlim(tt([1 end]));xlim(tt([1 end]));ylim([0 36]);xlim([tt_r(1) 60]);
    set(gca,'Xtick',0:20:120);
end

%% Main summary plot

%Summarize asymptote adaptation for all window sizes
[x_i,x_e,x_o] = deal(NaN(3,length(window),length(z.x.delay)/length(unique(z.x.delay))));

for n = 1:length(window)
    x.i{n}.values = nanmean(xii(setdiff(window{n},i_break),:),1);
    x.e{n}.values = nanmean(xee(setdiff(window{n},i_break),:),1);
    x.o{n}.values = nanmean(xoo(setdiff(window{n},i_break),:),1);
    x.i{n}.allvalues = xii(setdiff(window{n},i_break),:);
    x.e{n}.allvalues = xee(setdiff(window{n},i_break),:);
    x.o{n}.allvalues = xoo(setdiff(window{n},i_break),:);
    for m = 1:3
        x_i(m,n,:) = x.i{n}.values(z.x.delay==dList(m));
        x_e(m,n,:) = x.e{n}.values(z.x.delay==dList(m));
        x_o(m,n,:) = x.o{n}.values(z.x.delay==dList(m));
    end
end
%x_o = x_i+x_e;

for n = 1:3
    x.i{n}.values_by_latency = squeeze(x_i(:,n,:));
    x.e{n}.values_by_latency = squeeze(x_e(:,n,:));
    x.o{n}.values_by_latency = squeeze(x_o(:,n,:));
end

for n = 1:length(window)
    figure(102+n); hold on;
    set(102+n,'Position',[200,200,800,200])
    subplot(1,5,1:3); hold on;

    plot(dList2+2,nanmean(x_i(:,n,:),3),'linewidth',1,'color',[0.5 0.5 0.5]);
    plot(dList2-2,nanmean(x_e(:,n,:),3),'linewidth',1,'color',[0.5 0.5 0.5]);
    plot(dList2,nanmean(x_o(:,n,:),3),'linewidth',1,'color',[0.5 0.5 0.5]);

    for m = 1:3
        plot_force_errorbars_dot(dList2(m)+2,squeeze(x_i(m,n,:))',1,1,'color',c.blues(m,:));
        plot_force_errorbars_dot(dList2(m)-2,squeeze(x_e(m,n,:))',1,1,'color',c.reds(m,:));
        plot_force_errorbars_dot(dList2(m),squeeze(x_o(m,n,:))',1,1,'color',c.purps(m,:));
    end
    set(gca,'Xtick',dList2);set(gca,'Xticklabel',{'25ms','85ms','300ms'});
    xlim([0 325]);ylim([0 32]);

    [x.i{n}.p,x.e{n}.p,x.o{n}.p] = deal(NaN(3));
    for m = 1:3
        for q = (m+1):3
            [~,x.i{n}.p(m,q),~,x.i{n}.stats(m,q)] = ttest2(x_i(m,n,:),x_i(q,n,:),'tail','right');
            [~,x.e{n}.p(m,q),~,x.e{n}.stats(m,q)] = ttest2(x_e(m,n,:),x_e(q,n,:),'tail','left');
            [~,x.o{n}.p(m,q),~,x.o{n}.stats(m,q)] = ttest2(x_o(m,n,:),x_o(q,n,:),'tail','right');
        end
    end

    subplot(1,6,6); hold on;

    x_i_ = squeeze(x_i(:,n,:));

    diff_mean = [(mean(x_i_(1,:),'omitnan')-mean(x_i_(2,:),'omitnan')) (mean(x_i_(2,:),'omitnan')-mean(x_i_(3,:),'omitnan'))];
    diff_sem = [sqrt(nanvar(x_i_(1,:))/sum(~isnan(x_i_(1,:))) + nanvar(x_i_(2,:))/sum(~isnan(x_i_(2,:))))...
        sqrt(nanvar(x_i_(2,:))/sum(~isnan(x_i_(2,:))) + nanvar(x_i_(3,:))/sum(~isnan(x_i_(3,:))))];

    sens_mean = diff_mean./[6 21.5];
    sens_sem = diff_sem./[6 21.5];

    plot([1 2; 1 2],[sens_mean - sens_sem; sens_mean + sens_sem],'k');
    bar([1 2], sens_mean,0.75,'facecolor',c.blues(2,:));
    ylim([-0.2 1.8]);

    %Using F-test to compare sensitiities

    delays_ordered = [25*ones(1,34); 85*ones(1,34); 300*ones(1,34)]/10;
    [~,~,r_D] = regress(x_i_(:),[delays_ordered(:) delays_ordered(:)*0+1]);

    SSD_same = nansum(r_D.^2);
    SSD_separate = nansum(nansum(detrend(x_i_',0,'omitnan').^2));

    n_ = sum(~isnan(r_D));
    df_same = n_-2;
    df_separate = n_-3;

    F_i = (SSD_same-SSD_separate)/SSD_separate * df_separate/(df_same-df_separate);

    p_F_i = 1-fcdf(F_i,df_same-df_separate,df_separate);
    x.i{n}.sens_mean = sens_mean;
    x.i{n}.sens_sem = sens_sem;
    x.i{n}.diff_mean = diff_mean;
    x.i{n}.diff_sem = diff_sem;
    x.i{n}.F_i = F_i;
    x.i{n}.p_F_i = p_F_i;
    x.i{n}.n_F_i = n_;
end