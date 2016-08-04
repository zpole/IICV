function [im, scale] = readImage2(imagePath)

if ischar(imagePath)
  try
    im = imread(imagePath) ;
  catch
    error('Corrupted image %s', imagePath) ;
  end
else
  im = imagePath ;
end

scale = 1 ;
if (size(im,1) > 480)
  scale = 480 / size(im,1) ;
  im = imresize(im, scale) ;
  im = min(max(im,0),1) ;
end