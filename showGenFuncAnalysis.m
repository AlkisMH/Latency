function r = showGenFuncAnalysis(z_b,z_a,nIQR,c)

mvtDirs = -135:15:135;
dList=[0 85 300];   %latencies as coded in the experiment (experiment code would subtract 25 from the last two cases to calculate the extra delay)
dList2=[25 85 300]; %actual latency values corresponding to the above, for use in plotting
z_ab = combineLearningAndBaselineData(z_b,z_a);

%% Flipping and outlier rejection
z_a_flipped = (z_a.x.dat.g.imp).*shiftdim(sign(z_a.x.vmr),-1);
z_a_flipped_outlier_rej = NaN*z_a_flipped;

z_a_flipped_base1 = (z_a.x.dat.g.imp_base1).*shiftdim(sign(z_a.x.vmr),-1);
z_a_flipped_outlier_rej_base1 = NaN*z_a_flipped_base1;

z_a_exp_flipped = (z_a.x.dat.g.exp).*shiftdim(sign(z_a.x.vmr),-1);
z_a_exp_flipped_outlier_rej = NaN*z_a_exp_flipped;

z_a_exp_flipped_base1 = (z_a.x.dat.g.exp_base1).*shiftdim(sign(z_a.x.vmr),-1);
z_a_exp_flipped_outlier_rej_base1 = NaN*z_a_exp_flipped_base1;

z_a_tot_flipped = z_a_flipped+z_a_exp_flipped;
z_a_tot_flipped_outlier_rej = NaN*z_a_tot_flipped;

for m = 1:19
    for k = 1:3 %for each condition
        i_mk = z_a.x.delay==dList(k);
        z_a_flipped_outlier_rej_base1(:,m,i_mk) = removeOutliers(squeeze(z_a_flipped_base1(:,m,i_mk)),nIQR);
        z_a_exp_flipped_outlier_rej_base1(:,m,i_mk) = removeOutliers(squeeze(z_a_exp_flipped_base1(:,m,i_mk)),nIQR);

        z_a_flipped_outlier_rej(:,m,i_mk) = removeOutliers(squeeze(z_a_flipped(:,m,i_mk)),nIQR);
        z_a_exp_flipped_outlier_rej(:,m,i_mk) = removeOutliers(squeeze(z_a_exp_flipped(:,m,i_mk)),nIQR);
        z_a_tot_flipped_outlier_rej(:,m,i_mk) = removeOutliers(squeeze(z_a_tot_flipped(:,m,i_mk)),nIQR);
    end
end

z_b_flipped = (z_b.x.dat.g.imp).*shiftdim(sign(z_b.x.vmr),-1);
z_b_flipped_outlier_rej = NaN*z_b_flipped;

z_b_flipped_base1 = (z_b.x.dat.g.imp_base1).*shiftdim(sign(z_b.x.vmr),-1);
z_b_flipped_outlier_rej_base1 = NaN*z_b_flipped_base1;

for m = 1:19
    for k = 1:3 %for each condition
        i_mk = z_b.x.delay==dList(k);
        z_b_flipped_outlier_rej_base1(:,m,i_mk) = removeOutliers(squeeze(z_b_flipped_base1(:,m,i_mk)),nIQR);
        z_b_flipped_outlier_rej(:,m,i_mk) = removeOutliers(squeeze(z_b_flipped(:,m,i_mk)),nIQR);
    end
end

%% Calculate %outliers rejected

clear nq_b_before nq_b_after nq_a_before nq_a_after;
for i = 1:z_b.nsub
    i_ = z_b.x.participant==i;
    q_ = [z_b_flipped_base1(:,:,i_); z_b_flipped(:,:,i_)]; nq_b_before(i) = sum(isnan(q_(:)));
    q_ = [z_b_flipped_outlier_rej_base1(:,:,i_); z_b_flipped_outlier_rej(:,:,i_)]; nq_b_after(i) = sum(isnan(q_(:)));
end
for i = 1:z_a.nsub
    i_ = z_a.x.participant==i;
    q_ = [z_a_flipped_base1(:,:,i_); z_a_flipped(:,:,i_)]; nq_a_before(i) = sum(isnan(q_(:)));
    q_ = [z_a_flipped_outlier_rej_base1(:,:,i_); z_a_flipped_outlier_rej(:,:,i_)]; nq_a_after(i) = sum(isnan(q_(:)));
end

%re-write, but numel for Exp. 1b is 684 and for Exp. 1a is 456
%with base1:Exp. 1b:855; Exp. 1a: 570;
show_sterr_p(100*[nq_b_before/855 nq_a_before/570]);
show_sterr_p(100*[nq_b_after/855 nq_a_after/570]);
r.nq_b_before = nq_b_before;
r.nq_a_before = nq_a_before;
r.nq_b_after = nq_b_after;
r.nq_a_after = nq_a_after;

%% Average across trials
impG_a=squeeze(nanmean(z_a_flipped_outlier_rej,1));
impG_b=squeeze(nanmean(z_b_flipped_outlier_rej,1));
impG_ab = [impG_b, impG_a];

expG_a=squeeze(nanmean(z_a_exp_flipped_outlier_rej,1));

%Baseline
impB1G_a=squeeze(nanmean(z_a_flipped_outlier_rej_base1,1));
impB1G_b=squeeze(nanmean(z_b_flipped_outlier_rej_base1,1));
impB1G_ab = [impB1G_b, impB1G_a];

expB1G_a=squeeze(nanmean(z_a_exp_flipped_outlier_rej_base1,1));

%% Estimate local and global components

xx_a_expl = expG_a-expB1G_a;
%This won't be really used anywhere, the fitting parameters are
%recalculated below based on implicit from both Exp. 1a and 1b plus
%explicit from Exp. 1a

b_expl_a_param = NaN(2,size(xx_a_expl,2));
for k = 1:size(xx_a_expl,2)
    b_expl_a_param(:,k) = nlinfit(mvtDirs',xx_a_expl(:,k),@gaussian_fit_fixed,[12 4]);
end

%Implicit generalization
xx_ab = impG_ab-impB1G_ab;

clear b_impl_ab_param b_impl_exp_ab_shifted_param b_impl_exp_ab_shifted_param2;

bell_shape = gaussian_fit_fixed([1,0],mvtDirs);
b_impl_ab_param = NaN(2,size(xx_ab,2));
err = NaN(19,size(xx_ab,2));
for i = 1:size(xx_ab,2)
    %simplify: just regress onto the shapes
    [b_impl_ab_param(:,i),~,err(:,i)] = regress(xx_ab(:,i),[bell_shape(:),bell_shape(:)*0+1]);
    rsq_.rsq_0(i) = 1 - var(err(:,i))/var(xx_ab(:,i));
end

%compare across conditions
[r.i.exp_ab.local.p,r.i.exp_ab.global.p,r.e.exp_a.local.p,r.e.exp_a.global.p] = deal(NaN(3));
for m = 1:3
    for n = (m+1):3
        [~,r.i.exp_ab.local.p(m,n),~,r.i.exp_ab.local.stats(m,n)] = ttest2(b_impl_ab_param(1,z_ab.x.delay==dList(m)),b_impl_ab_param(1,z_ab.x.delay==dList(n)),'tail','right');
        [~,r.i.exp_ab.global.p(m,n),~,r.i.exp_ab.global.stats(m,n)] = ttest2(b_impl_ab_param(2,z_ab.x.delay==dList(m)),b_impl_ab_param(2,z_ab.x.delay==dList(n)),'tail','right');
        [~,r.e.exp_a.local.p(m,n),~,r.e.exp_a.local.stats(m,n)] = ttest2(b_expl_a_param(1,z_a.x.delay==dList(m)),b_expl_a_param(1,z_a.x.delay==dList(n)),'tail','left');
        [~,r.e.exp_a.global.p(m,n),~,r.e.exp_a.global.stats(m,n)] = ttest2(b_expl_a_param(2,z_a.x.delay==dList(m)),b_expl_a_param(2,z_a.x.delay==dList(n)),'tail','left');
    end
end

for m = 1:3
    r.i.exp_ab.local.values(m,:) = b_impl_ab_param(1,z_ab.x.delay==dList(m));
    r.i.exp_ab.global.values(m,:) = b_impl_ab_param(2,z_ab.x.delay==dList(m));
    r.e.exp_a.local.values(m,:) = b_expl_a_param(1,z_a.x.delay==dList(m));
    r.e.exp_a.global.values(m,:) = b_expl_a_param(2,z_a.x.delay==dList(m));
end

%% Main figure
figure(200);
set(200,'Position',[0 0 800 250])

xx_impl = impG_ab-impB1G_ab;
xx_expl = expG_a-expB1G_a; %though the baseline is essentiallly zero anyway

xoff = [0 -2 2];
for k = 1:3
    subplot(1,3,3); hold on;

    plot(-135:135,gaussian_fit_fixed(nanmean(b_expl_a_param(:,z_a.x.delay==dList(k)),2),-135:135),'color',0.5+0.5*c.reds(k,:),'linewidth',2);
    plot_force_errorbars_dot(mvtDirs+xoff(k),xx_expl(:,z_a.x.delay==dList(k)),1,1,'markersize',14,'linewidth',1,'markersize',14,'color',c.reds(k,:)); hold on;
    xlim([-150 150]);ylim([-1 32]);
    set(gca,'Xtick',-135:45:135);
    grid on;

    subplot(1,3,1); hold on;
    avg_params_tot = nanmean(b_expl_a_param(:,z_a.x.delay==dList(k)),2) + nanmean(b_impl_ab_param(:,z_ab.x.delay==dList(k)),2);
    plot(-135:135,gaussian_fit_fixed(avg_params_tot,-135:135),'color',0.5+0.5*c.purps(k,:),'linewidth',2);

    y_pred = gaussian_fit_fixed(avg_params_tot,-135:15:135);


    y_avg = nanmean(xx_impl(:,z_ab.x.delay==dList(k)),2) + nanmean(xx_expl(:,z_a.x.delay==dList(k)),2);
    y_sem = sqrt(nanvar(xx_impl(:,z_ab.x.delay==dList(k)),'',2)/34 + nanvar(xx_expl(:,z_a.x.delay==dList(k)),'',2)/16);

    plot([mvtDirs+xoff(k);mvtDirs+xoff(k)],[y_avg+y_sem y_avg-y_sem]','linewidth',1,'color',c.purps(k,:));
    plot(mvtDirs+xoff(k),y_avg,'.','markersize',14,'color',c.purps(k,:));
    xlim([-150 150]);ylim([-1 28]);
    set(gca,'Xtick',-135:45:135);
    grid on;

    r.o.Rsq(k) = 1-var(y_avg(:)-y_pred(:))/var(y_avg(:));

    subplot(1,3,2); hold on;
    plot(-135:135,gaussian_fit_fixed(nanmean(b_impl_ab_param(:,z_ab.x.delay==dList(k)),2),-135:135),'color',0.5+0.5*c.blues(k,:),'linewidth',2);
    plot_force_errorbars_dot(mvtDirs+xoff(k),xx_impl(:,z_ab.x.delay==dList(k)),1,1,'markersize',14,'linewidth',1,'markersize',14,'color',c.blues(k,:)); hold on;
    xlim([-150 150]);ylim([-1 32]);
    set(gca,'Xtick',-135:45:135);
    grid on;
end

%% Plot comparison of local/global components
figure(201);
set(201,'Position',[0 0 800 250])

subplot(1,3,1); hold on;
[b_tot_avg,b_tot_sem,b_impl_avg,b_expl_avg] = deal(NaN(2,3));
for m = 1:3
    b_tot_avg(1,m) = nanmean(b_impl_ab_param(1,z_ab.x.delay==dList(m)),2) + nanmean(b_expl_a_param(1,z_a.x.delay==dList(m)),2);
    b_tot_avg(2,m) = nanmean(b_impl_ab_param(2,z_ab.x.delay==dList(m)),2) + nanmean(b_expl_a_param(2,z_a.x.delay==dList(m)),2);

    b_tot_sem(1,m) = sqrt ( nanvar(b_impl_ab_param(1,z_ab.x.delay==dList(m)),'',2)/34 + nanvar(b_expl_a_param(1,z_a.x.delay==dList(m)),'',2)/16 );
    b_tot_sem(2,m) = sqrt ( nanvar(b_impl_ab_param(2,z_ab.x.delay==dList(m)),'',2)/34 + nanvar(b_expl_a_param(2,z_a.x.delay==dList(m)),'',2)/16 );
end
plot(dList2-1,b_tot_avg(1,:),'color',[0.5 0.5 0.5]);
plot(dList2+1,b_tot_avg(2,:),'color',[0.5 0.5 0.5]);
for m = 1:3
    plot(dList2(m)*[1 1]-1,b_tot_avg(1,m)+[b_tot_sem(1,m) -b_tot_sem(1,m)],'linewidth',1,'color',c.purps(m,:));
    plot(dList2(m)-1,b_tot_avg(1,m),'.','markersize',14,'color',c.purps(m,:));

    plot(dList2(m)*[1 1]+1,b_tot_avg(2,m)+[b_tot_sem(2,m) -b_tot_sem(2,m)],'linewidth',1,'color',c.purps(m,:));
    plot(dList2(m)+1,b_tot_avg(2,m),'.','markersize',14,'color',c.purps(m,:));

end

set(gca,'Xtick',dList2);set(gca,'Xticklabel',{'25','85','300'});
xlabel('Latency'); ylabel('Total local (degrees)');

xlim([0 325]); ylim([-2 25]);
grid on;

subplot(1,3,2); hold on;

for m = 1:3
    b_impl_avg(1,m) = nanmean(b_impl_ab_param(1,z_ab.x.delay==dList(m)),2);
    b_impl_avg(2,m) = nanmean(b_impl_ab_param(2,z_ab.x.delay==dList(m)),2);
end

plot(dList2,b_impl_avg(1,:),'color',[0.5 0.5 0.5]);
plot(dList2,b_impl_avg(2,:),'color',[0.5 0.5 0.5]);

for m = 1:3
    x_ = squeeze(b_impl_ab_param(1,z_ab.x.delay==dList(m))); x_ = x_(:)';
    plot_force_errorbars_dot(dList2(m),x_,1,1,'markersize',14,'color',c.blues(m,:));

    x_ = squeeze(b_impl_ab_param(2,z_ab.x.delay==dList(m))); x_ = x_(:)';
    plot_force_errorbars_dot(dList2(m),x_,1,1,'markersize',14,'color',c.blues(m,:));
end

set(gca,'Xtick',dList2);set(gca,'Xticklabel',{'25','85','300'});
xlabel('Latency'); ylabel('Implicit local (degrees)');

xlim([0 325]); ylim([-2 25]);
grid on;

subplot(1,3,3); hold on;
for m = 1:3
    b_expl_avg(1,m) = nanmean(b_expl_a_param(1,z_a.x.delay==dList(m)),2);
    b_expl_avg(2,m) = nanmean(b_expl_a_param(2,z_a.x.delay==dList(m)),2);
end

plot(dList2,b_expl_avg(1,:),'color',[0.5 0.5 0.5]);
plot(dList2,b_expl_avg(2,:),'color',[0.5 0.5 0.5]);

for m = 1:3
    x_ = squeeze(b_expl_a_param(1,z_a.x.delay==dList(m))); x_ = x_(:)';
    plot_force_errorbars_dot(dList2(m),x_,1,1,'markersize',14,'color',c.reds(m,:));

    x_ = squeeze(b_expl_a_param(2,z_a.x.delay==dList(m))); x_ = x_(:)';
    plot_force_errorbars_dot(dList2(m),x_,1,1,'markersize',14,'color',c.reds(m,:));
end

set(gca,'Xtick',dList2);set(gca,'Xticklabel',{'25','85','300'});
xlabel('Latency'); ylabel('Implicit local (degrees)');

xlim([0 325]); ylim([-2 25]);
grid on;

end