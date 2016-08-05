function filtrate_sift_by_slic_demo()

im = imread(fullfile(vl_root,'data','roofs1.jpg')) ;
im = im2single(im) ;
im = im(1:375,end-500+1:end,:) ;
  %slic
  segments = vl_slic(im, 30, 1, 'verbose') ;
  [sx,sy]=vl_grad(double(segments), 'type', 'forward') ;
  s = find(sx | sy) ;     %save the index of edge pixels
  %slic end
  
  h = size(im,1) ;
  features = getDenseSIFT(im) ;     %getdensesift（default）
  feature=filtrate_sift_by_slic(h,s,features, 'dsift');       %filtrate sift descriptor by slic edge
save testme.mat;