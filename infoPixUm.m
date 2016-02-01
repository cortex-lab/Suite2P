function um2pix=infoPixUm(totPixels,zoomMicro,microID)

% this is Sylvias version
% infoPixUm computes the number of pixels per um along the horizontal
% direction (xPU) vertical direction (yPU) approximately any direction
% (rPU)
%
% The inputs are:

% totPixels (number of pixels along the horizontal or vertical directions
% in the pre-registered movie)
% zoomMicro (zoom of the microscope)
% microID 'b'=b-scope, 'm'=MOM
%
% 2015.12.16 Mario Dipoppa - Created
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if zoomMicro==3 && microID=='b'
    xPU=totPixels/370;
    yPU=totPixels/340;
elseif zoomMicro==2 && microID=='b' 
    xPU=totPixels/520;
    yPU=totPixels/500;
elseif zoomMicro==4 && microID=='m'
    xPU=totPixels/117;
    yPU=totPixels/117;
elseif zoomMicro==3 && microID=='m'
    xPU=totPixels/155;
    yPU=totPixels/155;
end

rPU=sqrt(yPU^2+xPU^2)/sqrt(2);

um2pix.xPU=xPU;
um2pix.yPU=yPU;
um2pix.rPU=rPU;

end
