%Try to make my life easier when I import figures to the Illustrator
set(groot,'defaultAxesXColor',[0 0 0]);
set(groot,'defaultAxesYColor',[0 0 0]);
set(groot,'defaultAxesFontName','Arial');

% load crunched data
if ~exist('z_a'), z_a=load('exp1a_data.mat'); end
if ~exist('z_b'), z_b=load('exp1b_data.mat'); end

% reoganize feature data from separate cell arrays for each participant to a unified matrix for each feature
% And include only outward mvts in the reoganized data
z_a = reorganizeData(z_a,1);
z_b = reorganizeData(z_b,1);

% IQR-based outlier removal across participants
nIQR = 3;  % to be used later

ii_a = buildEpochIndexVectors(2);
ii_b = buildEpochIndexVectors(3);

%Move this to generalization function function
for i = 1:z_a.nsub
    z_a.tgt_file_all{i} = [];
    for k = 1:length(z_a.final_data{i}.tgt_file)
        z_a.tgt_file_all{i} = [z_a.tgt_file_all{i};z_a.final_data{i}.tgt_file{k}(2:end-1,:)];
    end
    if length(z_a.final_data{i}.tgt_file)==31
        z_a.tgt_file_all{i} = [z_a.tgt_file_all{i};NaN(60,20)];
    end
end

for i = 1:z_b.nsub
    z_b.tgt_file_all{i} = [];
    for k = 1:length(z_b.final_data{i}.tgt_file)
        z_b.tgt_file_all{i} = [z_b.tgt_file_all{i};z_b.final_data{i}.tgt_file{k}(2:end-1,:)];
    end
end

z_a.x.rawdat.aim_dir_gen = NaN*z_a.x.rawdat.aim_dir;
z_a.x.rawdat.aim_dir_gen0 = NaN*z_a.x.rawdat.aim_dir;

for k = 1:z_a.nsub
    mask = z_a.tgt_file_all{k}(1:2:end,5)==0;
    z_a.x.rawdat.aim_dir_gen(:,k) = z_a.x.rawdat.aim_dir(:,k).*mask./mask;
    z_a.x.rawdat.aim_dir_gen0(:,k) = z_a.x.rawdat.aim_dir0(:,k).*mask./mask;
end

z_b.x.rawdat.aim_dir_gen = NaN*z_b.x.rawdat.aim_dir;
z_b.x.rawdat.aim_dir_gen0 = NaN*z_b.x.rawdat.aim_dir0;

za_0 = organizeLearningAndBaselineData(z_a,ii_a,1);
zb_0 = organizeLearningAndBaselineData(z_b,ii_b,1);

za_0 = organizeGeneralizationData(za_0,ii_a,1);
zb_0 = organizeGeneralizationData(zb_0,ii_b,1);

z_0 = combineLearningAndBaselineData(zb_0,za_0);

%% Color scheme for plots
blue1 = [ 0,0,205]./255;
blue2 = [ 30,144,255]./255;
blue3 = [64,224,208]./255;
red1 = [139,0,0]./255;
red2 = [255,0,0]./255;
red3 = [240,128,128]./255;
purp1 = [75,0,130]./255;
purp2 = [128,0,128]./255;
purp3 = [153,50,204]./255;

c.blues = [blue1;blue2;blue3];
c.reds = [red1;red2;red3];
c.purps = [purp1;purp2;purp3];

%% Learning analysis
analysis_window = (101:120)+16;
analysis_window_relearn = (41:60)+136;
analysis_window_combined = [analysis_window analysis_window_relearn];

x_ab_0 = plotLearning(z_0, nIQR, c, {analysis_window,analysis_window_relearn,analysis_window_combined});

%% Numbers on paper
clc;
n = 1; %1: learning; 2: relearning; 3: combined
disp('------------------');
disp('-----Overall------');
disp('------------------');
show_sterr_p(x_ab_0.o{n}.values_by_latency);
x_ab_0.o{n}.p
[x_ab_0.o{n}.stats(1,2).tstat,x_ab_0.o{n}.stats(1,2).df]
[x_ab_0.o{n}.stats(1,3).tstat,x_ab_0.o{n}.stats(1,2).df]
disp('------------------');
disp('-----Implicit-----');
disp('------------------');
show_sterr_p(x_ab_0.i{n}.values_by_latency);
x_ab_0.i{n}.p
[x_ab_0.i{n}.stats(1,2).tstat,x_ab_0.i{n}.stats(1,2).df]
[x_ab_0.i{n}.stats(1,3).tstat,x_ab_0.i{n}.stats(1,3).df]
disp('------------------');
disp('-----Explicit-----');
disp('------------------');
show_sterr_p(x_ab_0.e{n}.values_by_latency);
x_ab_0.e{n}.p
[x_ab_0.e{n}.stats(1,2).tstat,x_ab_0.e{n}.stats(1,2).df]
[x_ab_0.e{n}.stats(1,3).tstat,x_ab_0.e{n}.stats(1,3).df]

%% Baseline comparison (without flipping)

%Note: here, I'm to compare baselines based on the sign of the PRECEDING
%superset. But, since the preceding sign is always the opposite of the
%current sign, I can select based on current sign

nanmean(x_ab_0.i{1}.base_signed(z_0.x.order>1 & sign(z_0.x.vmr)>0)) -...
    nanmean(x_ab_0.i{1}.base_signed(z_0.x.order>1 & sign(z_0.x.vmr)<0))

[~,p,~,stats] = ttest2(x_ab_0.i{1}.base_signed(z_0.x.order>1 & sign(z_0.x.vmr)>0),x_ab_0.i{1}.base_signed(z_0.x.order>1 & sign(z_0.x.vmr)<0))

nanmean(x_ab_0.e{1}.base_signed(z_0.x.order>1 & sign(z_0.x.vmr)>0)) -...
    nanmean(x_ab_0.e{1}.base_signed(z_0.x.order>1 & sign(z_0.x.vmr)<0))

[~,p,~,stats] = ttest2(x_ab_0.e{1}.base_signed(z_0.x.order>1 & sign(z_0.x.vmr)>0),x_ab_0.e{1}.base_signed(z_0.x.order>1 & sign(z_0.x.vmr)<0))

nanmean(x_ab_0.o{1}.base_signed(z_0.x.order>1 & sign(z_0.x.vmr)>0)) -...
    nanmean(x_ab_0.o{1}.base_signed(z_0.x.order>1 & sign(z_0.x.vmr)<0))

[~,p,~,stats] = ttest2(x_ab_0.o{1}.base_signed(z_0.x.order>1 & sign(z_0.x.vmr)>0),x_ab_0.o{1}.base_signed(z_0.x.order>1 & sign(z_0.x.vmr)<0))

%% Generalization analysis
r_0 = showGenFuncAnalysis(zb_0,za_0,nIQR,c);

%% Check outlier rejection - % rejected etc.

yyy = [z_0.x.dat.head_dir_peak_sp0_base;z_0.x.dat.head_dir_peak_sp0];

clear z_bef z_aft;
for i = 1:z_0.nsub
    i_ = z_0.x.participant==i;
    q_ = isnan(yyy(:,i_));
    z_bef(i) = sum(q_(:));
    z_aft(i) = sum(x_ab_0.NNr(i_));
end

%Combine with generalization data
total_trials_examined = [3*196+855*ones(1,18) 2*196+570*ones(1,24)];
total_rej_before = z_bef + [r_0.nq_b_before r_0.nq_a_before];
total_rej_after = z_aft + [r_0.nq_b_after r_0.nq_a_after];

show_sterr_p(100*total_rej_before./total_trials_examined);
show_sterr_p(100*total_rej_after./total_trials_examined);
show_sterr_p(100*total_rej_after./total_trials_examined-100*total_rej_before./total_trials_examined);

%% Check movement/acceleration times

ii_test_a = [horzcat(ii_a.Base1{:}),horzcat(ii_a.Learn1{:}),horzcat(ii_a.genAll{:}),horzcat(ii_a.genBase1{:})];

za_0.x.rawdat.acc_time = NaN(1345,za_0.nsub);
for k = 1:za_0.nsub
    q_ = 5*(za_0.final_data{k}.peak_sp_P(1:2:end) - za_0.final_data{k}.i_start_recalc(1:2:end));
    za_0.x.rawdat.acc_time(1:length(q_),k) = q_;
end

ii_test_b = [horzcat(ii_b.Base1{:}),horzcat(ii_b.Learn1{:}),horzcat(ii_b.genAll{:}),horzcat(ii_b.genBase1{:})];
zb_0.x.rawdat.acc_time = NaN(1856,zb_0.nsub);
for k = 1:zb_0.nsub
    zb_0.x.rawdat.acc_time(:,k) = 5*(zb_0.final_data{k}.peak_sp_P(1:2:end) - zb_0.final_data{k}.i_start_recalc(1:2:end));
end

acc_time = [nanmean(zb_0.x.rawdat.acc_time(ii_test_b,:)) nanmean(za_0.x.rawdat.acc_time(ii_test_a,:))];

show_sterr_p(acc_time);
