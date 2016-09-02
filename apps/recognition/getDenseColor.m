function features = getDenseColor(im)

%options.scale                     = 1;
%%options.sigma_scale               = 0.6;
%opts.deltax                       = 120;
%opts.deltay                       = 90;
%options.color                     = 1;
%options.nbins                     = 3;
%options.patchsize                 = 2;
%options.norm                      = 1;
%%options.clamp                     = 0.2;
[descr, frame] = denseCOLOR(im);  %[descr, frame] = denseCOLOR(im,options);
frame([1,2],:)=frame([2,1],:);
frame(7,:)=frame(7,:)/2;
frame(3:6,:)=[];
descr = single(descr);
features.descr = descr;
features.frame = frame;