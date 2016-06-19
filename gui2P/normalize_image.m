function mimg = normalize_image(mimg)

mimg = mimg - my_conv2(mimg, 5, [1 2]);

mimg = mimg ./ my_conv2(abs(mimg), 10, [1 2]);
