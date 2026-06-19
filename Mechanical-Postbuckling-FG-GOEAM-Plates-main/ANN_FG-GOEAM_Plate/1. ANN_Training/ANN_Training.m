%% ========================================================================
%  Author:        Vu Ngoc Viet Hoang
%  Affiliation:   Bach Khoa Phu Tho Technology Co., Ltd
%  Location:      Ho Chi Minh City, Viet Nam
%  Email:         vnvhoang2610@gmail.com
%  Date:          March 2025
% =========================================================================

clear all;
clc;

% Load data from Excel file
data = readtable('Uniaxial_X.xlsx');

inputs = data{:, 1:8}'; 
targets = data{:, 9:11}'; 

%% Bayesian Regularization 
% hiddenLayerSizes = [15 15]; 
% net = fitnet(hiddenLayerSizes, 'trainbr'); 

%% Levenberg-Marquardt 
hiddenLayerSizes = [15 15]; 
net = fitnet(hiddenLayerSizes, 'trainlm');

%% Scaled conjugate gradient
% hiddenLayerSizes = [15 15]; 
% net = fitnet(hiddenLayerSizes, 'trainscg'); 


% Set transfer functions for hidden layers
for i = 1:length(hiddenLayerSizes)
    net.layers{i}.transferFcn = 'tansig'; % Hyperbolic tangent sigmoid transfer function tansig (alternatively, 'logsig')
end

% Set transfer function for the output layer
net.layers{end}.transferFcn = 'purelin'; 

%% Specify how data is divided into training, testing, and validation sets
net.divideFcn = 'dividerand';
net.divideMode = 'sample'; 

%% Set training ratios for Bayesian Regularization- NOTE: "trainbr" does not require a validation dataset
% net.divideParam.trainRatio = 70/100; 
% net.divideParam.testRatio = 30/100; 

%% Set training ratios for Levenberg-Marquardt and Scaled conjugate gradient
net.divideParam.trainRatio = 70/100; 
net.divideParam.valRatio = 15/100; 
net.divideParam.testRatio = 15/100; 

%% Set training parameters
net.trainParam.epochs = 2000; 
net.trainParam.goal = 1e-7; 

%% Train the neural network
[net, tr] = train(net, inputs, targets);

%% Save the trained network and training record to a .mat file
save('LM_Uniaxial_X.mat', 'net', 'tr');

