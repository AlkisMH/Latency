function z = organizeLearningAndBaselineData(z,ii,zero_flag)
% excerpt out the data from VMR learing and baseline epochs
Ns = length(z.final_data);
Nb = length(ii.Learn1);
if zero_flag
    fn = {'head_dir_peak_sp0','aim_dir0','peak_sp','pert','vfb','tgt_angle'};
else
    fn = {'head_dir_peak_sp','aim_dir','peak_sp','pert','vfb','tgt_angle'};
end
for s=1:Ns,
    for b=1:Nb,
        k2 = (s-1)*Nb + b;
        for k=1:length(fn),
            zz = z.x.rawdat.(fn{k})(ii.Learn1{b},s); Nz(k,k2)=length(zz);
            z.x.dat.(fn{k})(1:Nz(k,k2),k2) = zz;
            if (k==1)|(k==2),
                zz = z.x.rawdat.(fn{k})(ii.Base1{b},s); Nz(k,k2)=length(zz);
                z.x.dat.([fn{k},'_base'])(1:Nz(k,k2),k2) = zz;
            end
        end
        z.x.order(k2)       = b;
        z.x.participant(k2) = s;
        z.x.delay(k2)       = median( z.final_data{s}.tgt_file{12*b}(:,2) );
        z.x.vmr(k2)         = median( z.x.dat.pert(:,k2) );
        z.x.trainTgtDir(k2) = round(nanmean(z.x.dat.tgt_angle(:,k2))/15)*15;
        if isnan(z.x.trainTgtDir(k2)), keyboard; end
    end
end
z.NzBlock = Nz;
z.x.dat.tot = z.x.dat.(fn{1});
z.x.dat.exp = z.x.dat.(fn{2});
z.x.dat.imp = z.x.dat.tot - z.x.dat.(fn{2});
z.x.dat.tot_base = z.x.dat.([fn{1},'_base']);
z.x.dat.exp_base = z.x.dat.([fn{2},'_base']);
z.x.dat.imp_base = z.x.dat.tot_base - z.x.dat.([fn{2},'_base']);
end