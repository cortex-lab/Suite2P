function TiffWriter(image,fname,bitspersamp)
% written by Adam Packer (I think).

t = Tiff(fname,'w');
tagstruct.ImageLength = size(image,1);
tagstruct.ImageWidth = size(image,2);
tagstruct.Photometric = Tiff.Photometric.MinIsBlack;
if bitspersamp==16
    tagstruct.BitsPerSample = 16;
end
if bitspersamp==32
    tagstruct.BitsPerSample = 32;
end
tagstruct.SampleFormat = Tiff.SampleFormat.Int;
tagstruct.SamplesPerPixel = 1;
tagstruct.RowsPerStrip = 256;
tagstruct.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;
tagstruct.Software = 'MATLAB';
t.setTag(tagstruct);
t.write(image(:,:,1));
for i=2:size(image,3)
    t.writeDirectory();
    t.setTag(tagstruct);
    t.write(image(:,:,i));
end
t.close();