im = imread(fullfile(vl_root,'data','roofs1.jpg')) ;
%im = im2single(im) ;
%im = im(1:128,end-128+1:end,:) ;
%im = im2single(im) ;
%im = vl_imsmooth(im, 0) ;

[a,b] = denseCOLOR(im);

save testdcolor.mat