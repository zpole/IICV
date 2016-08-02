im = imread(fullfile(vl_root,'data','roofs1.jpg')) ;
im2 = im2single(im) ;
%im = im2single(im) ;
%im = im(1:128,end-128+1:end,:) ;
%im = im2single(im) ;
%im = vl_imsmooth(im, 0) ;

%slic
segments = vl_slic(im2, 30, 1, 'verbose') ;
[sx,sy]=vl_grad(double(segments), 'type', 'forward') ;
s = find(sx | sy) ;     %save the index of edge pixels
%slic end
  
h = size(im2,1) ;

[features.dcolor , features.infodcolor] = denseCOLOR(im);
feature=filtrate_sift_by_slic(h,s,features, 'dcolor'); 
save testdcolor.mat