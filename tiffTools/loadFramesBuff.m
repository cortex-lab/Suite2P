function [frames, headers] = loadFramesBuff(tiff, firstIdx, lastIdx, stride, temp_file)
%loadFrames Loads the frames of a Tiff file into an array (Y,X,T)
%   MOVIE = loadFrames(TIFF, [FIRST], [LAST], [STRIDE], []) loads
%   frames from the Tiff file specified by TIFF, which should be a filename
%   or an already open Tiff object. Optionallly FIRST, LAST and STRIDE
%   specify the range of frame indices to load.

    if nargin>4 && ~isempty(temp_file)
        if ~isequal(tiff, temp_file) % do not copy if already copied
            % in case copying fails (server hangs)
            iscopied = 0;
            firstfail = 1;
            while ~iscopied
                try
                    copyfile(tiff,temp_file);
                    iscopied = 1;
                    if ~firstfail
                        fprintf('  succeeded!\n');
                    end
                catch
                    if firstfail
                        fprintf('copy tiff failed, retrying...');
                    end
                    firstfail = 0;
                    pause(10);
                end
            end
            tiff = temp_file;
        end
        info = imfinfo(temp_file); % get info after copying
        if isnan(lastIdx)
            lastIdx = length(info); % get number of frames
        end
    end

    % initChars = overfprintf(0, 'Loading TIFF frame ');
    warning('off', 'MATLAB:imagesci:tiffmexutils:libtiffWarning');

    warningsBackOn = onCleanup(@() warning('on', 'MATLAB:imagesci:tiffmexutils:libtiffWarning'));

    if ischar(tiff) || isstring(tiff)
        tiff = Tiff(tiff, 'r');
        closeTiff = onCleanup(@() close(tiff));
    end

    if nargin < 2 || isempty(firstIdx)
        firstIdx = 1;
    end

    if nargin < 3 || isempty(lastIdx)
        lastIdx = nFramesTiff(tiff);
    end

    if nargin < 4 || isempty(stride)
        stride = 1;
    end

    if nargout > 1
        loadHeaders = true;
    else
        loadHeaders = false;
    end

    w = tiff.getTag('ImageWidth');
    h = tiff.getTag('ImageLength');
    dataClass = class(read(tiff));
    nFrames = ceil((lastIdx - firstIdx + 1)/stride);
    frames = zeros(h, w, nFrames, dataClass);
    if loadHeaders
        headers = cell(1, nFrames);
    end

    setDirectory(tiff, firstIdx);
    for t = 1:nFrames
        frames(:, :, t) = read(tiff);

        if loadHeaders
            headers{t} = getTag(tiff, 'ImageDescription');
            try
                headers{t} = [headers{t} getTag(tiff,'Software')];
            catch
            end
        end

        if t < nFrames
            for i = 1:stride
                nextDirectory(tiff);
            end
        end
    end
end
