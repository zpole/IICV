im = imread(fullfile(vl_root,'data','roofs1.jpg')) ;
im2 = im2single(im) ;

%slic
%segments = vl_slic(im2, 30, 1, 'verbose') ;
%[sx,sy]=vl_grad(double(segments), 'type', 'forward') ;
%s = find(sx | sy) ;     %save the index of edge pixels
%slic end
  
%h = size(im2,1) ;

feature = getDenseColor(im);
%feature=filtrate_sift_by_slic(h,s,features, 'dcolor'); 
save testdcolor.mat