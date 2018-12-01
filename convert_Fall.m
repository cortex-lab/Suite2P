% python to matlab suite2p
% python version saves Fall file, this converts Fall to a format compatible
% with the matlab GUI

dat = load('Fall.mat');

%%
Fcell = {dat.F};
FcellNeu = {dat.Fneu};
ops = dat.ops;
ops.Ly = single(dat.ops.Ly);
ops.Lx = single(dat.ops.Lx);

ops.mimg1 = ops.meanImg;
if isfield(ops, 'meanImg_chan2')
	ops.mimgRED = ops.meanImg_chan2;
end
sp = {dat.spks};
clear stat
flds = fieldnames(dat.stat{1});
for n = 1:length(dat.stat)
	for j = 1:length(flds)
		stat(n).(flds{j}) = dat.stat{n}.(flds{j});
	end
	stat(n).ipix = int64(ops.Ly)*(stat(n).xpix) + stat(n).ypix + 1;
	stat(n).ipix = stat(n).ipix(:);
	stat(n).mimgProjAbs = 0;
	stat(n).cmpct = stat(n).compact;
	stat(n).aspect_ratio = double(stat(n).aspect_ratio);
	%stat(n) = dat.stat{n};
end
ops.yrange = [1:dat.ops.Ly];
ops.xrange = [1:dat.ops.Lx];
ops.yrange_crop = dat.ops.yrange(1)+1:dat.ops.yrange(end);
ops.xrange_crop = dat.ops.xrange(1)+1:dat.ops.xrange(end);

ops.Vcorr = zeros(dat.ops.Ly, dat.ops.Lx);
ops.Vcorr(ops.yrange_crop, ops.xrange_crop) = dat.ops.Vcorr;


%%


save('F.mat','-v7.3','Fcell','FcellNeu','ops','sp','stat');
