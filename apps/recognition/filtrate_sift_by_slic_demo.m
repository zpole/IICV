function filtrate_sift_by_slic_demo()

im = imread(fullfile(vl_root,'data','cat.jpeg')) ;
im = im2single(im) ;
%im = im(1:375,end-500+1:end,:) ;
  %slic
  segments = vl_slic(im, 50, 0.1, 'verbose') ;
  [sx,sy]=vl_grad(double(segments), 'type', 'forward') ;
  s = find(sx | sy) ;     %save the index of edge pixels
  %slic end
  
  h = size(im,1) ;
  features = getDenseSIFT(im) ;     %getdensesift???default???
  
 
  sel = vl_colsubset(1:size(features.descr,2), 5120) ;
  feature.frame = features.frame(:,sel) ;
  feature.descr = features.descr(:,sel) ;
  feature.contrast = features.contrast(:,sel) ;
  
  feature=filtrate_sift_by_slic(h,s,features, 'dsift');       %filtrate sift descriptor by slic edge
  
  figure(3) ; clf ;
  imp = im ;
  imp([s s+numel(im(:,:,1)) s+2*numel(im(:,:,1))]) = 0 ;
  imagesc(imp);axis image off ; hold on ;

   h1 = vl_plotframe(feature.frame) ;
  set(h1,'color','y','linewidth',1) ;
save testme.mat;