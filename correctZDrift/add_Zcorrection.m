
function add_Zcorrection(ops)

%
ops.Zexpt                   = getOr(ops, 'Zexpt', ops.expts(end)+1);
ops.Zplanes                 = getOr(ops, 'Zplanes', 200);
ops.Zchannels               = getOr(ops, 'Zchannels', 1);
ops.Zzoom                   = getOr(ops, 'Zzoom',     1);
ops.ZstackSavePath          = 'D:\DATA\zStacks';

% registration options
ops.NiterPrealign       = getOr(ops, 'NiterPrealign', 5);
ops.fracImgPreAlign     = getOr(ops, 'fracImgPreAlign', 1);

% check if Z-stack already registered
fname = sprintf('stack_%s_%s.mat', ops.mouse_name, ops.date);
fs    = dir(fullfile(ops.ZstackSavePath, fname));

% register Z-stack
if isempty(fs)
   registerZstack(ops); 
else
    disp('z-stack already registered');
end

%% align mean planes to z-stack

[PtoZ, MimgZ, Zpatch, Bimg] = affinePlanestoZ(ops);
ipl = [];
for j = 1:length(PtoZ)
    if ~isempty(PtoZ{j})
        ipl = [ipl j];
    end
end

%% align frames across recording to z-stack, compute position zpos
% write plane in z-stack at that zpos to binary file
zpos = writeZtoBin(ops, PtoZ, MimgZ);

%% load binary files and extract z-stack signals
disp('extracting signals from binary files...');
for iplane = ipl
    planefile = sprintf('%s/F_%s_%s_plane%d.mat', ops.ResultsSavePath, ...
        ops.mouse_name, ops.date, iplane);
    dat = load(planefile);
    opsZ = ops;
    opsZ.RegFile              = sprintf('%s_Z.bin',dat.ops.RegFile(1:end-4));
    [~, ~, FcellZ, FcellNeuZ] = extractSignalsSurroundNeuropil(opsZ, dat.stat);   
    
    % save FcellZ and FcellNeuZ
    dat.FcellZ    = FcellZ;
    dat.FcellNeuZ = FcellNeuZ;
    
    save(planefile, '-struct', 'dat')
end

