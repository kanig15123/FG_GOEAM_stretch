%% ========================================================================
%  Author:        Vu Ngoc Viet Hoang
%  Affiliation:   Bach Khoa Phu Tho Technology Co., Ltd
%  Location:      Ho Chi Minh City, Viet Nam
%  Email:         vnvhoang2610@gmail.com
%  Date:          March 2025
% =========================================================================

%% ------------------------- Input Parameters ----------------------------
var1 = 15;
var2 = 1.5; 
var3 = 350;  
var4 = 85;  
var5 = 2;   
var6 = 5; 
var7 = 10;
var8 = 15; 

% Load the trained ANN model
load('BR_Uniaxial_U.mat', 'net', 'tr');

% Define the new input data
newInput = [var1; var2; var3; var4; var5; var6; var7; var8];

predictedOutput = net(newInput);

% Assign output values to coefficients
Lambda0 = predictedOutput(1);
Lambda1 = predictedOutput(2);
Lambda2 = predictedOutput(3);

%% ----------------------------- Plotting --------------------------------
num_points = 100; 
Wn = linspace(0, 1, num_points); % Normalized W displacement

Lambda_n = Lambda0 + Lambda1 .* Wn + Lambda2 .* (Wn.^2);

% Create the plot
figure;
plot(Wn, Lambda_n, 'b-', 'LineWidth', 3.5);

grid on; 
box on;
set(gca, 'LineWidth', 2, 'FontSize', 26, 'TickLabelInterpreter', 'latex');
pbaspect([4 3 1]);

xlabel('$\overline{W}$', 'FontSize', 28, 'Interpreter', 'latex');  
ylabel('$\overline{\lambda}$', 'FontSize', 28, 'Interpreter', 'latex');  
title('\textbf{\boldmath$\mathbf{U\!-\!W}_{\mathrm{Gr}}$ Plate Under Uiaxial Compression}', ...
      'FontSize', 30, 'FontWeight', 'bold', 'Interpreter', 'latex');
