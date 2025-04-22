function z = reorganizeData(z,zero_flag)
% reoganize feature data from separate cell arrays for each participant to a unified matrix for each feature
% And include only outward mvts in the reoganized data

if zero_flag
    fn = {'head_dir_peak_sp','aim_dir','head_dir_peak_sp0','aim_dir0','peak_sp','pert','vfb','tgt_angle'};
else
    fn = {'head_dir_peak_sp','aim_dir','peak_sp','pert','vfb','tgt_angle'};
end
for k=1:length(fn)
    for s=1:length(z.final_data)
        zz = z.final_data{s}.(fn{k})(1:2:end); Nz(k,s)=length(zz);
        z.x.rawdat.(fn{k})(1:Nz(k,s),s) = zz;
        %Note that, if shorter than maximum, the end will be appended with
        %zeros - initialize with NaNs? (we don't really use the last block
        %though)
    end
end

% Also add vector variables - here I'm only doing velocity for now
for s=1:length(z.final_data)
    [all_vels,all_tgt_vels,all_tgt_pos] = deal(NaN(121,length(z.final_data{s}.Vel_tab)));
    for i = 1:length(z.final_data{s}.Vel_tab)
        istart = z.final_data{s}.i_start_recalc(i);
        if isfinite(istart)
            v_ = [NaN(20,1);z.final_data{s}.Vel_tab{i};NaN(100,1)];
            all_vels(:,i) = v_(istart:(istart+120))*2.54/96;

            %Note that tgt-direction referenced movements are already
            %normalized
            v_ = [NaN(20,1);z.final_data{s}.Vel_tgt{i};NaN(100,1)];
            all_tgt_vels(:,i) = v_(istart:(istart+120));

            v_ = [NaN(20,1);z.final_data{s}.Pos_tgt{i};NaN(100,1)];
            all_tgt_pos(:,i) = v_(istart:(istart+120));
        else
            %otherwise nothing gets written, so we have NaNs as initialized
        end
    end
    zz = all_vels(:,1:2:end);
    z.x.rawdat.all_vels(:,1:length(zz),s) = zz;

    zz = all_tgt_vels(:,1:2:end);
    z.x.rawdat.all_tgt_vels(:,1:length(zz),s) = zz;

    zz = all_tgt_pos(:,1:2:end);
    z.x.rawdat.all_tgt_pos(:,1:length(zz),s) = zz;
end

%% Add extra info: block lengths, ITI, acceleration time, peak speed
%Initializations
[z.x.rawdat.ITI,z.x.rawdat.acc_time,z.x.rawdat.mvt_time,z.x.rawdat.max_vel] = deal(NaN*z.x.rawdat.(fn{1}));

%Block lengths
z.block_lengths = (cellfun(@length,z.final_data{1}.tgt_file)-2)/2;

%ITI
for s=1:length(z.final_data)
    qqq = cell2mat(z.final_data{s}.start_and_end_times);
    aa1 = qqq(1:2:end,1);
    i_ = [find(diff(aa1)<0); length(aa1)];
    i1_ = [0; i_];
    trial_start_times = NaN*aa1;
    for i = 1:length(i_)-1
        trial_start_times((i1_(i)+1):i1_(i+1)) = ...
            aa1((i1_(i)+1):i1_(i+1))-aa1(i1_(i+1))+z.final_data{s}.t_last_trial(i);
    end
    z.x.rawdat.ITI(1:length(trial_start_times),s) = [0;diff(trial_start_times)];
end

for s=1:length(z.final_data)

    %keyboard;
    %Acceleration time
    q_ = 5*(z.final_data{s}.peak_sp_P(1:2:end) - z.final_data{s}.i_start_recalc(1:2:end));
    z.x.rawdat.acc_time(1:length(q_),s) = q_;

    %Movement times (NEEDS DOUBLECHECKING!!)
    %mvt_times_raw = vertcat(z.final_data{s}.start_and_end_times{:,2})-vertcat(z.final_data{s}.start_and_end_times{:,1});
    mvt_times_raw = cellfun(@length,z.final_data{s}.Vel_tab)*5;
    q_ = mvt_times_raw(1:2:end);
    z.x.rawdat.mvt_time(1:length(q_),s) = q_;

    % Maximum velocity
    max_vels = NaN*mvt_times_raw;
    for i = 1:length(max_vels)
        i_vmax = z.final_data{s}.peak_sp_P(i);
        if isfinite(i_vmax)
            max_vels(i) = z.final_data{s}.Vel_tab{i}(z.final_data{s}.peak_sp_P(i))*2.54/96; %note conversion from pixels/s to cm/s
        end
    end
    q_ = max_vels(1:2:end);
    z.x.rawdat.max_vel(1:length(q_),s) = q_;
end
z.NzWhole=Nz;
end