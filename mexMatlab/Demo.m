clear, clc;

%% Compiling and Choosing input files(image domain in [0,255], uint8)
mex -O SCFmex.cpp;
image = imread('..\testimg\test.jpg');

%% Setting the parameters
k=2000;
sigmar = 0.2;
iteration = 5;
step= 2;

%% Filtering using soft clustering filtering
[smoothing(:,:,1), smoothing(:,:,2), smoothing(:,:,3)] = SCFmex(image, k, sigmar, iteration, step);

%% Showing result
figure, imshow(uint8(smoothing));