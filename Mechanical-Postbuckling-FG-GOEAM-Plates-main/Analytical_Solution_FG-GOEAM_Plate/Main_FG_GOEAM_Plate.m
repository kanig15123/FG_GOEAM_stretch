%% ========================================================================
%  Author:        Vu Ngoc Viet Hoang
%  Affiliation:   Bach Khoa Phu Tho Technology Co., Ltd
%  Location:      Ho Chi Minh City, Viet Nam
%  Email:         vnvhoang2610@gmail.com
%  Date:          March 2025
% =========================================================================

clear all
clc

%% ------------------------- Input Parameters ----------------------------
var1 = 15;
var2 = 1.5; 
var3 = 350;  
var4 = 85;  
var5 = 2;   
var6 = 5; 
var7 = 10;
var8 = 15; 

%% Uniaxial Compression
Chi1=-1;
Chi2=0;

%% Biaxial Compression
% Chi1=-1;
% Chi2=-1;

%% Material and geometrical parameters 
h=0.005;     
a=var1*h;
b=var2*a;
T=var3;
HGr=var4/100;
WGr=var5/100; 

% Mode numbers
m=1;
n=1;

% Constants for HSDT formulation
c1 = 4/h^2/3;
c2 = 3*c1;

%% ------------------ Material Properties for FG-GOEAM -------------------
FG_GOEAM; 

%% ------------------------- Kerr Foundation -----------------------------
Kun=var6;
Ksn=var7;
Kln=var8;

Ku=Kun*Dm/(a^4);
Kl=Kln*Dm/(a^4);
Ks=Ksn*Dm/(a^2);

K1=Kl*Ku/(Kl+Ku);
K2=Ks*Ku/(Kl+Ku);

%% ------------------- Analytical Solution Parameters --------------------
Parameters_Analytical_Solution;

%% Eqs. (38-39)
Lambda0= lamda0 .*a.^2./(ECu.*h.^3);
Lambda1= lamda1 .*a.^2./(ECu.*h.^3);
Lambda2= lamda2 .*a.^2./(ECu.*h.^3);

%% ----------------------------- Plotting --------------------------------

num_points = 100; 
Wn = linspace(0, 1, num_points); % Normalized W displacement

Lambda_n = Lambda0 + Lambda1 .* Wn + Lambda2 .* (Wn.^2);

% Create the plot
figure;
plot(Wn, Lambda_n, 'k-', 'LineWidth', 3.5);

grid on; 
box on;
set(gca, 'LineWidth', 2, 'FontSize', 26, 'TickLabelInterpreter', 'latex');
pbaspect([4 3 1]);

xlabel('$\overline{W}$', 'FontSize', 28, 'Interpreter', 'latex');  
ylabel('$\overline{\lambda}$', 'FontSize', 28, 'Interpreter', 'latex');  
title('\textbf{\boldmath$\mathbf{U\!-\!W}_{\mathrm{Gr}}$ Plate Under Uiaxial Compression}', ...
      'FontSize', 30, 'FontWeight', 'bold', 'Interpreter', 'latex');


