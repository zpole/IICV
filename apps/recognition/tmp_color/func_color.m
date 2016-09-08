function [features_color,color_info] = func_color(im)
%function features_color = func_color(im)
%im = imread(fullfile(vl_root,'data','roofs1.jpg')) ;

 options.scale                     = 1;
 %options.sigma_scale               = 0.6;
 opts.deltax                       = 120;
 opts.deltay                       = 90;
 options.color                     = 1;
 options.nbins                     = 3;
 options.patchsize                 = 2;
 options.norm                      = 1;
 %options.clamp                     = 0.2;
 
 [dcolor , infodcolor]             = denseCOLOR(im, options ); 
 
 save col.mat
 
 features_color = dcolor;
 color_info = infodcolor;