% Find the number of frames in a tiff file
%
% This function reports number of frames in a tiff file by
% jumping over seekInterval frames until the end of the file.
% This is in general faster than Matlab's imfinfo function.
%
%  USAGE
%   n = nFrames(tiff, seekInterval)
%   tiff            Path to tiff file or a handler of an already opened file.
%   seekInterval    Optional number of frames to jump over.
%                   Default is 1000.
%   n               Number of frames
%
function n = nFrames(tiff, seekInterval)
    if nargin < 2
        seekInterval = 1000;
    end
    %keep guessing until we seek too far
    guess = seekInterval;
    overSeeked = false;
    n = 0;

    if ischar(tiff) || isstring(tiff)
        tiff = Tiff(tiff, 'r');
        closeTiff = onCleanup(@() close(tiff));
    end

    while ~overSeeked
      try
        tiff.setDirectory(guess);
        guess = 2*guess; %double the guess
      catch ex
        overSeeked = true; %we tried to seek past the last directory
      end
    end
    %when overseeking occurs, the current directory/frame will be the last one
    n = tiff.currentDirectory;
end

