function [nPlanes, nChannels] = nPlanesChannels(fname)

[~, header] = loadFramesBuff(fname, 1, 1, 1);

hh=header{1};
str = hh(strfind(hh, 'channelsSave = '):end);
ind = strfind(str, 'scanimage');
ch = str2num(str(16 : ind(1)-1));
nChannels = length(ch);

fastZEnable = sscanf(hh(findstr(hh, 'fastZEnable = '):end), 'fastZEnable = %d');
fastZDiscardFlybackFrames = sscanf(hh(findstr(hh, 'fastZDiscardFlybackFrames = '):end), 'fastZDiscardFlybackFrames = %d');
if isempty(fastZDiscardFlybackFrames)
    fastZDiscardFlybackFrames = 0;
end
stackNumSlices = sscanf(hh(findstr(hh, 'stackNumSlices = '):end), 'stackNumSlices = %d');

nPlanes = 1;
if fastZEnable
    nPlanes = stackNumSlices+fastZDiscardFlybackFrames;
end