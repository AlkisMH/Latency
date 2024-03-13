function z = organizeGeneralizationData(z,ii,zero_flag)
% excerpt out the data from VMR learing and baseline epochs
% instead of a 120x1 LC for each block, make a 6x19 GF for each
Ns = length(z.final_data);
Nb = length(ii.Learn1);
tgtList = [-135:15:135]; Nt = length(tgtList);
if zero_flag
    fn = {'head_dir_peak_sp','head_dir_peak_sp0','aim_dir','aim_dir_gen','aim_dir0','aim_dir_gen0','pert','vfb','tgt_angle'};
else
    fn = {'head_dir_peak_sp','aim_dir','aim_dir_gen','pert','vfb','tgt_angle'};
end
% rmfield(z.x.dat,"g");
for s=1:Ns, %subjects
    for b=1:Nb, %different conditions
        k2 = (s-1)*Nb + b;
        for k=1:length(fn), %metrics
            i = ii.genAll{b};
            zt = round(z.x.rawdat.tgt_angle(i,s)) - z.x.trainTgtDir(k2);
            zt = round(360*(zt/360 - round(zt/360)));
            zz = z.x.rawdat.(fn{k})(i,s);
            for kk=1:Nt, %directions
                [kk,tgtList(kk)];
                zz1 = [zz(zt==tgtList(kk));NaN];
                Nz(k2,kk)=length(zz1);
                %try
                z.x.dat.g.(fn{k})(1:12,kk,k2) = zz1(1:12);
                %catch
                %    keyboard;
                %end
            end
            if (zero_flag & k<7)||k<4, %only compute baseline for some of the metrics (baseline not meaningful for all metrics!)
                i = ii.genBase1{b}; %use preceding baseline for each condition
                zt = round(z.x.rawdat.tgt_angle(i,s)) - z.x.trainTgtDir(k2);
                zt = round(360*(zt/360 - round(zt/360)));
                zz = z.x.rawdat.(fn{k})(i,s);
                for kk=1:Nt,
                    zzt = zz(zt==tgtList(kk));
                    try
                    z.x.dat.g.([fn{k},'_base1'])(1:length(zzt),kk,k2) = zzt;
                    catch
                    keyboard;
                    end
                end
                i = ii.genBase2{b}; %use baselines from all conditions from the same participant
                zt = round(z.x.rawdat.tgt_angle(i,s)) - z.x.trainTgtDir(k2);
                zt = round(360*(zt/360 - round(zt/360)));
                zz = z.x.rawdat.(fn{k})(i,s);
                for kk=1:Nt,
                    zzt = zz(zt==tgtList(kk)); Nz2=length(zzt); zzt(Nz2+1:length(i)/19) = NaN;
                    z.x.dat.g.([fn{k},'_base2'])(:,kk,k2) = zzt;
                end
            end
        end
    end
end
% Nz-12

if zero_flag
    z.x.dat.g.imp = z.x.dat.g.head_dir_peak_sp0 - z.x.dat.g.aim_dir0;
    z.x.dat.g.exp = z.x.dat.g.aim_dir_gen0;
    z.x.dat.g.imp_base1 = z.x.dat.g.head_dir_peak_sp0_base1 - z.x.dat.g.aim_dir0_base1;
    z.x.dat.g.exp_base1 = z.x.dat.g.aim_dir0_base1;
    z.x.dat.g.imp_base2 = z.x.dat.g.head_dir_peak_sp0_base2 - z.x.dat.g.aim_dir0_base2;
    z.x.dat.g.exp_base2 = z.x.dat.g.aim_dir0_base2;
else
    z.x.dat.g.imp = z.x.dat.g.head_dir_peak_sp - z.x.dat.g.aim_dir;
    z.x.dat.g.exp = z.x.dat.g.aim_dir_gen;
    z.x.dat.g.imp_base1 = z.x.dat.g.head_dir_peak_sp_base1 - z.x.dat.g.aim_dir_base1;
    z.x.dat.g.exp_base1 = z.x.dat.g.aim_dir_base1;
    z.x.dat.g.imp_base2 = z.x.dat.g.head_dir_peak_sp_base2 - z.x.dat.g.aim_dir_base2;
    z.x.dat.g.exp_base2 = z.x.dat.g.aim_dir_base2;
end
end