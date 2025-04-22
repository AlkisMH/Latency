function ii = buildEpochIndexVectors(nBlocks)
% Note this works for the two Experiments (1a and 1b), because they share 
% the same duration familiarization and training epochs
% Initial submission version 3/23/2024
% 1/31/2025: I expanded in order to save indices for washout data as well
o = cumsum([380,511,511]);
for k =1:nBlocks,
    ii.Base1{k}   = o(k) + [1:16]; % [381:396];
    ii.Learn1{k}  = o(k) + [16+[1:120], 16+120+114+[1:60]]; % [[397:516], [631:690]];
    ii.genAll{k}  = o(k) + [16+120+[1:114], 16+120+114+60+[1:114]];
    ii.genBase1{k} = o(k) + [-56:0]; %generalization baseline 1: separate for each block
    ii.genBase2{k} = o(1:nBlocks) + [-56:0]'; ii.genBase2{k} = ii.genBase2{k}(:)'; %generalization baseline 2: same for all two or three blocks
    ii.Wash1{k} = ii.genAll{k}(end) + [1:30];
end
end