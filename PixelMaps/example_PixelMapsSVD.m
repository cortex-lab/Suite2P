addpath(genpath('D:\CODE\GitHub\Suite2P\'))

iplane = 3;
mname = 'M150610_MP020';
date_exp = '2015-11-10';
block = 5;

%% this file includes the starts of my stimuli, but you might have them in another form
allinfo_path = sprintf('//zserver\\Data\\expInfo\\%s\\%s\\%d\\%s_%d_%s_allinfo.mat', mname, date_exp, block, ...
    date_exp, block, mname);
load(allinfo_path)

%% this collects event types into a cell array 'ton'
ievt = allinfo.eventMatrix.parsedEvents.eventsStimAndRepeat(:,1);
ft1 = allinfo.ft1{iplane};

evtuniq = unique(ievt);
ton = cell(length(evtuniq), 1);
for i = 1:length(evtuniq)
   ton{i} = ft1(ievt==evtuniq(i));
end

%% path to example SVD file is on zserver 
% (hint: it's faster to copy this file locally first)

root = '\\zserver\Lab\Share\Marius\Examples';
% root = sprintf('D:\\DATA\\F\\%s\\%s', mname, date_exp); % this is my normal local path
fname = sprintf('SVD_%s_%s_plane%d.mat', mname, date_exp, iplane);
SVDfilepath = fullfile(root, fname);

ops.ExpFiltTau  = 30;
ops.ntf         = 60; % timepoints for the filter
ops.method      = 'sta'; % deconv or sta
ops.nSVDforPixelMaps = 250;

[R, ops] = PixelMapsSVD(ton, block, SVDfilepath, ops);

%% find maximum of responses (this block will be replaced by von Mises fits)
[amps, bestphi]  = max(reshape(R, [], size(R,3)), [], 2);
amps             = reshape(amps, size(R,1), size(R,2));
bestphi          = reshape(bestphi, size(R,1), size(R,2));

%% plot HSV map of responses. Value is max response. 
[Ly, Lx, nstims] = size(R);

Hue = bestphi/size(R,3);
Sat = 1 * ones(Ly, Lx);

V   = .4 * (amps/mean(amps(:)));
V = min(1, max(0, V));
% V = .85 * ones(Ly,Lx);
% V(Sat<.2) = 1;

close all
figure('Position', [10, 10, 1400, 1400])
Ibig = hsv2rgb(cat(3, Hue, Sat, V));

% Ibig = imresize(Ibig, 1, 'lanczos3'); %resize to increase resolution if necessary

imagesc(Ibig)

axis square
set(gcf, 'Color', 'w')
axis off

% export_fig 'D:/FIG/NicePics/oripix.pdf' -q101

%%