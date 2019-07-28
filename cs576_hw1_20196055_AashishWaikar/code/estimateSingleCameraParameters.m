function [cameraParams] = estimateSingleCameraParameters(imagePoints, boardSize, patchSize, imageSize)
% This function will estimate camera parameters (intrinsic, extrinsic) from
% checkerboard image points.

% Zhang's method consists of 5 parts
% 1. Estimate homography from checkerboard plane to screen space.
% 2. Calculate B matrix by solving Vb = 0.
% 3. Extract intrinsic parameters from B matrix.
% 4. Calculate extrinsic parameters from intrinsic parameters and homography.
% 5. Refine parameters using the maximum likelihood estimation.

% Function inputs:
% - 'imagePoints': positions of checkerboard points in a screen space.
% - 'boardSize': the number of horizontal, vertical patchs in the checkerboard.
% - 'patchSize': the size of the checkerboard patch in mm.
% - 'imageSize': the size of the checkerboard image in pixels.

% Function outputs:
% - 'cameraParams': a camera parameter includes intrinsic and extrinsic.

numView = size(imagePoints, 3);
numVerticalPatch = boardSize(1) - 1;
numHorizontalPatch = boardSize(2) - 1;
numCorner = size(imagePoints, 1);

%% Estimate a homography (appendix A)
% Generate checkerboard world points
worldPoints = zeros(size(imagePoints,1), size(imagePoints,2));

% Fill worldPoints (positions of checkerboard corners)
% ----- Your code here (1) ----- (slide 6)

%for i=0:numCorner-1
%        worldPoints(i+1,1) = mod(i,9)*30;
%        worldPoints(i+1,2) = mod(i,6)*30;
%end
for i=1:numHorizontalPatch
    for j=1:numVerticalPatch
        worldPoints(numVerticalPatch*(i-1)+j,1)=(i-1)*patchSize;
        worldPoints(numVerticalPatch*(i-1)+j,2)=(j-1)*patchSize;
    end
end


% Build L matrix
L = zeros(2 * numCorner, 9, numView);
for i=1:numView
    j=1;
    while j<=numCorner
        L(2*(j-1)+1,1,i)=-worldPoints(j,1);
        L(2*(j-1)+1,2,i)=-worldPoints(j,2);
        L(2*(j-1)+1,3,i)=-1;
        L(2*(j-1)+1,7,i)=worldPoints(j,1)*imagePoints(j,1,i);
        L(2*(j-1)+1,8,i)=worldPoints(j,2)*imagePoints(j,1,i);
        L(2*(j-1)+1,9,i)=imagePoints(j,1,i);
        L(2*(j-1)+2,4,i)=-worldPoints(j,1);
        L(2*(j-1)+2,5,i)=-worldPoints(j,2);
        L(2*(j-1)+2,6,i)=-1;
        L(2*(j-1)+2,7,i)=worldPoints(j,1)*imagePoints(j,2,i);
        L(2*(j-1)+2,8,i)=worldPoints(j,2)*imagePoints(j,2,i);
        L(2*(j-1)+2,9,i)=imagePoints(j,2,i);
        j=j+1;
    end
end


% Fill L matrix
% ----- Your code here (2) ----- (slide 13)



% Calculate a homography using SVD
homography = zeros(3,3,numView);

% Fill homography matrix
% ----- Your code here (3) ----- (slide 15)
for v=1:numView
    [U,S,V]=svd(L(:,:,v));
    %s=svd(L(:,:,v));
    %S
    %d=size(s)
    %ssv=s(d(1));
    %c=1;
    %while S(c,c)~=ssv
    %    c=c+1;
    %end
    h=V(:,9);
    h=h/h(9);
    homography(1,1,v)=h(1);
    homography(1,2,v)=h(2);
    homography(1,3,v)=h(3);
    homography(2,1,v)=h(4);
    homography(2,2,v)=h(5);
    homography(2,3,v)=h(6);
    homography(3,1,v)=h(7);
    homography(3,2,v)=h(8);
    homography(3,3,v)=h(9);
end
%homography
%% Solve closed-form (section 3.1)
V = zeros(2 * numView, 6);
b = zeros(6, 1);

% Fill V matrix and calculate b vector
% ----- Your code here (4) ----- (slide 19, 23)
for i=1:numView
    %size(h(1,1,i))
    %h(1,2,i)
    V(2*i-1,1)=homography(1,1,i)*homography(1,2,i);
    V(2*i-1,2)=homography(1,1,i)*homography(2,2,i) + homography(2,1,i)*homography(1,2,i);
    V(2*i-1,3)=homography(1,1,i)*homography(3,2,i) + homography(3,1,i)*homography(1,2,i);
    V(2*i-1,4)=homography(2,1,i)*homography(2,2,i);
    V(2*i-1,5)=homography(2,1,i)*homography(3,2,i) + homography(3,1,i)*homography(2,2,i);
    V(2*i-1,6)=homography(3,1,i)*homography(3,2,i);
    
    V(2*i,1)=homography(1,1,i)*homography(1,1,i) - homography(1,2,i)*homography(1,2,i);
    V(2*i,2)=homography(1,1,i)*homography(2,1,i) + homography(2,1,i)*homography(1,1,i) - homography(1,2,i)*homography(2,2,i) + homography(2,2,i)*homography(1,2,i);
    V(2*i,3)=homography(1,1,i)*homography(3,1,i) + homography(3,1,i)*homography(1,1,i) - homography(1,2,i)*homography(3,2,i) + homography(3,2,i)*homography(1,2,i);
    V(2*i,4)=homography(2,1,i)*homography(2,1,i) - homography(2,2,i)*homography(2,2,i);
    V(2*i,5)=homography(2,1,i)*homography(3,1,i) + homography(3,1,i)*homography(2,1,i) - homography(2,2,i)*homography(3,2,i) + homography(3,2,i)*homography(2,2,i);
    V(2*i,6)=homography(3,1,i)*homography(3,1,i) - homography(3,2,i)*homography(3,2,i);
end
%V
[U1,S1,V1]=svd(V);
% S1
% s1=svd(V);
% d1=size(s1);
% ssv1=s1(d1);
% c1=1;
% while S1(c1,c1)~=ssv1
%     c1=c1+1;
% end

b = V1(:,6);

%% Extraction of the intrinsic parameters from matrix B (appendix B)

% ----- Your code here (5) ----- (slide 24)
v0 = (b(2)*b(3)-b(1)*b(5))/(b(1)*b(4)-b(2)*b(2));  % modify this line
lambda = b(6) - (b(3)*b(3)+v0*(b(2)*b(3)-b(1)*b(5)))/b(1);  % modify this line
alpha = sqrt(lambda/b(1));  % modify this line
beta = sqrt(lambda*b(1)/(b(1)*b(4)-b(2)*b(2)));  % modify this line
gamma = -b(2)*alpha*alpha*beta/lambda;  % modify this line
u0 = gamma*v0/beta - b(3)*alpha*alpha/lambda;  % modify this line

p=[v0,lambda,alpha,beta,gamma,u0];


%% Estimate initial RT (section 3.1)
Rt = zeros(3, 4, numView);
K_init = zeros(3,3);
K_init(1,1)=alpha;
K_init(1,2)=gamma;
K_init(1,3)=u0;
K_init(2,2)=beta;
K_init(2,3)=v0;
K_init(3,3)=1;
% Fill Rt matrix
% ----- Your code here (6) ----- (slide 25, 26)
for i=1:numView
    h1=homography(:,1,i);
    h2=homography(:,2,i);
    h3=homography(:,3,i);
    norm1=norm(inv(K_init)*h1);
    norm2=norm(inv(K_init)*h2);
    lambda1=(1/norm1 + 1/norm2)/2;
    r1=lambda1*inv(K_init)*h1;
    r2=lambda1*inv(K_init)*homography(:,2,i);
    r3=cross(r1,r2);
    t=lambda1*inv(K_init)*h3;
    R=zeros(3,3);
    R(:,1)=r1;
    R(:,2)=r2;
    R(:,3)=r3;
    
    [U,S,V] = svd(R);
    R = U*V';
    
    Rt(:,1,i)=R(:,1);
    Rt(:,2,i)=R(:,2);
    Rt(:,3,i)=R(:,3);
    Rt(:,4,i)=t;
end

%K_init
%Rt

%% Maximum likelihood estimation (section 3.2)
options = optimoptions(@lsqnonlin, 'Algorithm', 'levenberg-marquardt', ...
    'TolX', 1e-32, 'TolFun', 1e-32, 'MaxFunEvals', 1e64, ...
    'MaxIter', 1e64, 'UseParallel', true);

% Build initial x value as x0
% ----- Your code here (7) ----- (slide 29)


% 5 for intrinsic
% 3 for translation, 3 for rotation, total 6 for each checkerboard image
x0 = zeros(5 + 6 * size(imagePoints, 3), 1);  % modify this line
x0(1)=alpha;
x0(2)=beta;
x0(3)=gamma;
x0(4)=u0;
x0(5)=v0;

Rvec = zeros(3,numView);
for i=1:numView
    Rvec(:,i)=rotationMatrixToVector(Rt(:,1:3,i)');
end

for i=0:numView-1
    x0(5+6*i+1)=Rt(1,4,i+1);
    x0(5+6*i+2)=Rt(2,4,i+1);
    x0(5+6*i+3)=Rt(3,4,i+1);
    x0(5+6*i+4)=Rvec(1,i+1);
    x0(5+6*i+5)=Rvec(2,i+1);
    x0(5+6*i+6)=Rvec(3,i+1);
end

% Non-least square optimization
% Read [https://mathworks.com/help/optim/ug/lsqnonlin.html] for more information
[objective] = @(x) func_calibration(imagePoints, worldPoints, x);

[x_hat, ~, ~, ~, ~] = lsqnonlin(objective,x0,[],[],options);


%% Build camera parameters
rvecs = zeros(numView, 3);
tvecs = zeros(numView, 3);
K = [1, 0, 0
     0, 1, 0
     0, 0, 1];

% Extract intrinsic matrix K, rotation vectors and translation vectors from x_hat
% ----- Your code here (8) -----
K(1,1)=x_hat(1);
K(2,2)=x_hat(2);
K(1,2)=x_hat(3);
K(1,3)=x_hat(4);
K(2,3)=x_hat(5);

for i=1:numView
    tvecs(i,1)=x_hat(5+6*(i-1)+1);
    tvecs(i,2)=x_hat(5+6*(i-1)+2);
    tvecs(i,3)=x_hat(5+6*(i-1)+3);
    rvecs(i,1)=x_hat(5+6*(i-1)+4);
    rvecs(i,2)=x_hat(5+6*(i-1)+5);
    rvecs(i,3)=x_hat(5+6*(i-1)+6);
end

% Generate cameraParameters structure
cameraParams = cameraParameters('IntrinsicMatrix', K', ...
    'RotationVectors', rvecs, 'TranslationVectors', tvecs, ...
    'WorldPoints', worldPoints, 'WorldUnits', 'mm', ...
    'imageSize', imageSize) ; 


reprojected_errors = zeros(size(imagePoints));

% Uncomment this line after you implement this function to calculate
% reprojection errors of your camera parameters.
 reprojected_errors = imagePoints - cameraParams.ReprojectedPoints;

cameraParams = cameraParameters('IntrinsicMatrix', K', ...
    'RotationVectors', rvecs, 'TranslationVectors', tvecs, ...
    'WorldPoints', worldPoints, 'WorldUnits', 'mm', ...
    'imageSize', imageSize, 'ReprojectionErrors', reprojected_errors) ; 