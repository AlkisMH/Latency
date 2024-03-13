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
    end
end
z.NzWhole=Nz;
end