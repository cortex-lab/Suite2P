function n = nFramesTiff(tiff)
%nFrames find the number of frames in the Tiff

%keep guessing until we seek too far
guess = 2001;
overSeeked = false;

if ischar(tiff)
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

