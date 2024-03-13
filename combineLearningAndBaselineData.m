function z = combineLearningAndBaselineData(z1,z2)
% excerpt out the data from VMR learing and baseline epochs
fn = fieldnames(z1.x.dat);

for k=1:length(fn)
    z.x.dat.(fn{k}) = [z1.x.dat.(fn{k}), z2.x.dat.(fn{k})];
end
z.x.order = [z1.x.order, z2.x.order];
z.x.participant = [z1.x.participant, max(z1.x.participant) + z2.x.participant];
z.x.delay = [z1.x.delay, z2.x.delay];
z.x.vmr = [z1.x.vmr, z2.x.vmr];
z.x.trainTgtDir = [z1.x.trainTgtDir, z2.x.trainTgtDir];
z.NzBlock = [z1.NzBlock, z2.NzBlock];
z.nsub = length(unique(z.x.participant));
end
