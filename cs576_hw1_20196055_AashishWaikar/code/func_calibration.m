function [objective] = func_calibration(imagePoints, worldPoints, x)
% Objective function to minimize eq.10 in Zhang's paper. 
% Size of input variable x is 5+6*n where n is number of checkerboard 
% images. An intrinsic matrix can be reconstructed from first five
% parameters, and the extrinsic matrix can be reconstructed from remain
% parameters.

% You should fill the variable hat_m which contains reprojected positions 
% of checkerboard points in screen coordinate.

% Function inputs:
% - 'imagePoints': positions of checkerboard points in a screen space.
% - 'worldPoints': positions of checkerboard points in a model space.
% - 'x': parameters to be optimized.

% Function outputs:
% - 'objective': difference of estimated values and real values.
    
numView = size(imagePoints,3);
numCor = size(imagePoints,1);
hat_m = zeros(size(imagePoints));

K_in = zeros(3,3);
K_in(1,1)=x(1);
K_in(2,2)=x(2);
K_in(1,2)=x(3);
K_in(1,3)=x(4);
K_in(2,3)=x(5);
K_in(3,3)=1;
for i=1:numView
    for j=1:numCor
        %Ho=zeros(3,3);
        rt=zeros(3,3);
        rot_mat = rotationVectorToMatrix(x(5+6*(i-1)+4:5+6*(i-1)+6))';
        rt(:,1:2)=rot_mat(:,1:2);
        rt(:,3)=x(5+6*(i-1)+1:5+6*(i-1)+3);
        H=K_in*rt;
        wh=zeros(3);
        wh(1:2)=worldPoints(j,:);
        wh(3)=1;
        ih=H*wh;
        hat_m(j,:,i)=ih(1:2)/ih(3);
    end
end
% ----- Your code here (9) -----


objective = imagePoints - hat_m;