function [depthMap, disparityMap] = estimateDepth(leftImage, rightImage, stereoParameters)
% This function estimate disparity and depth values from left and right
% images. You should calculate disparty map first and then convert the
% disparity map to depth map using left camera parameters.

% Function inputs:
% - 'leftImage': rectified left image.
% - 'rightImage': rectified right image.
% - 'stereoParameters': stereo camera parameters.

% Function outputs:
% - 'depth': depth map of left camera.
% - 'disparity': disparity map of left camera.

leftImageGray = rgb2gray(im2double(leftImage));
rightImageGray = rgb2gray(im2double(rightImage));

translation = stereoParameters.TranslationOfCamera2;
baseline = norm(translation);
focalLength = stereoParameters.CameraParameters1.FocalLength(1);

disparityMap = zeros(size(leftImageGray));
depthMap = zeros(size(leftImageGray));

%size(leftImage)
% ----- Your code here (10) -----

max_d = 583; 
min_d = 0; 
 
  window = ones(13);
  window = window ./ numel(window);
  %kernel

  I1 = leftImageGray;
  I2 = rightImageGray;

  d_vals = min_d : max_d;
  num_d = length(d_vals);
  C = zeros(size(I1,1), size(I1,2), num_d); %cost volume
  

  for d = 1 : length(d_vals);
    I2t = imtranslate(I2, [d 0]);
    C(:,:,d) = abs(I1 - I2t); 
    C(:,:,d) = imfilter(C(:,:,d), window);

  end

  [C_min, D] = min(C, [], 3);
  disparityMap = D + min_d;
  for i=1:size(disparityMap,1)
      for j=1:size(disparityMap,2)
          depthMap(i,j) = focalLength*baseline/disparityMap(i,j);
      end
  end

end