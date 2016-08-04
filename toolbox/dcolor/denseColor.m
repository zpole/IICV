% denseCOLOR compute histograms of color projection on a regular dense grid
% 
% Usage
% ------
% 
% [dcolor , infodcolor] = denseCOLOR(I , [options] );
% 
% 
% Inputs
% -------
% 
% I                                     Input image (ny x nx x [3]) in UINT8 format. 
% 
% options
% 	   scale                          Scaling vector (1 x nscale). Extract descriptors at different scaling of the image (default scale = [1]).
% 	   sigma_scale                    Scaling factor to obtain the standard deviation of the Gaussian filter (sigma = sigma_scale/scale)(default sigma_scale = 0.6)
% 	   deltax                         Division step in the x-axis for the grid (default deltax = floor(nx*min(scale))) 
% 	   deltay                         Division step in the y-axis for the grid (default deltay = floor(ny*min(scale)))
%        color                          0 : force gray-scale (dimcolor = 1, default), 1 : RGB (dimcolor = 3), 2 : nRGB (dimcolor = 3), 3 : Opponent (dimcolor = 3), 
%                                       4 : nOpponent (dimcolor = 2), 5 : Hue (dimcolor = 1)
% 	   nbins                          Number of bins for histograms (default nbins = 32)
% 	   patchsize                      Size of the patch where the descriptor is computed (default patchsize = nbins/2 )	  
% 	   norm                           Normalization : norm = 0 <=> no normalization, norm = 1 <=> v=v/(sum(v)+epsi), norm = 2 <=> v=v/sqrt(sum(v?+epsi?, 
% 	                                  norm = 3 <=> v=sqrt(v/(sum(v)+epsi)) , norm = 3 <=> L2-clamped (default norm = 1)
% 	   clamp                          Clamping value (default clamp = 0.2)
% 
% 
% Outputs
% -------
% 
% dcolor                                COLOR descriptors (nbins x nb_pts) where nb_pts = deltax*deltay*nscale*dimcolor
% infodcolor                            COLOR descriptors informations(7 x nb_pts)   where nb_pts = deltax*deltay*nscale*dimcolor
%                                       infodcolor(1,i) = y
% 									  infodcolor(2,i) = x
% 									  infodcolor(3,i) = scale
% 									  infodcolor(4,i) = color
% 									  infodcolor(5,i) = nyscale
% 									  infodcolor(6,i) = nxscale
% 									  infodcolor(7,i) = currentpatchsize
% 									  
% 
% 
% Example 1
% ---------
% 
% 
% clear all;
% I                                    = imread('E:/ImageAudio/Data/graz02/cars/carsgraz_001.bmp');
% % I = imresize(I);
% 
% options.scale                        = 1;
% % options.sigma_scale                  = 0.6;
% % opts.deltax = 100;
% % opts.deltay = 100;
% options.color                        = 1;
% options.nbins                        = 8;
% % options.patchsize                    = ;
% options.norm                         = 0;
% % options.clamp                        = 0.2;
% 
% [dcolor , infodcolor]                 = denseCOLOR(I , options ); 
% 
% figure(1)
% imagesc(dcolor)
% 
% figure(2)
% imagesc(I)
% colormap(gray)
% hold on
% plot(infodcolor(2 , :) , infodcolor(1 , :) , 'r+')
% hold off


% Example 2
% ---------
% 
% 
% I                                    = imread('image_0174.jpg');
% 
% 
% options.scale                        = [1];
% options.sigma_scale                  = 0.6;
% options.deltax                       = 10;
% options.deltay                       = 10;
% options.nbins                        = 16;
% options.patchsize                    = 16;
% options.norm                         = 1;
% 
% 
% [dcolor , infodcolor]                = denseCOLOR(I , options ); 
% 
% half                                 = options.patchsize/2;
% 
% figure(1)
% imagesc(dcolor)
% 
% xr                                   = [infodcolor(2, :)-half ; infodcolor(2, :)-half ; infodcolor(2, :)+ half ; infodcolor(2, :)+ half ; infodcolor(2, :)-half] + 1.5;
% yr                                   = [infodcolor(1, :)-half ; infodcolor(1, :)+half ; infodcolor(1, :)+ half ; infodcolor(1, :)- half ; infodcolor(1, :)-half] + 1.5;
% 
% 
% figure(2)
% imagesc(I)
% colormap(gray)
% hold on
% plot(infodcolor(2 , :)+1.5 , infodcolor(1 , :)+1.5 , 'r+')
% plot(xr , yr , 'b')
% hold off
% 
% 
% To compile
% ----------
% 
% mex  -output denseCOLOR.dll denseCOLOR.c
% 
% mex  -g -output denseCOLOR.dll denseCOLOR.c
% 
% mex  -f mexopts_intel10.bat -output denseCOLOR.dll denseCOLOR.c
% 
% 
% Author : Sébastien PARIS : sebastien.paris@lsis.org
% -------  Date : 02/04/2010
% 
% 
% References :  
% ---------        