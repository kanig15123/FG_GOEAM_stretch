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
h=0.002;
a=10*h;           
b=a;

% Mode numbers
m=1;
n=1;

% Constants for HSDT formulation
c1 = 4/h^2/3;
c2 = 3*c1;

% Power law indices of FGM
p = 1;  

%% ------------------ Material Properties for FGM -------------------------
FGM; 

%% ---------------------- Pasternak Foundation ----------------------------
Dm=Em*h^3/12/(1-Num^2);

Kwn=100; % Winkler spring
Gpn=10;  % Pasternak shear layer

Kw=Kwn*Dm/a^4;
Gp=Gpn*Dm/a^2;

K1=Kw;
K2=Gp;

%% Uniaxial Compression
% Chi1=-1;
% Chi2=0;

%% Biaxial Compression
Chi1=-1;
Chi2=-1;

%% ------------------- Analytical Solution Parameters --------------------
Parameters_Analytical_Solution;

%% Eq. (36) 
Lamda0=(-1).*(S2+V11+P11.*V14+P21.*V15).*(Chi1.*V18+Chi2.*V19).^(-1);

Lamda1=(-1).*h.*(V12+P12.*V14+P22.*V15+P11.*V16+P21.*V17).*(Chi1.*V18+Chi2.*V19).^(-1);

Lamda2=(-1).*h.^2.*(V13+P12.*V16+P22.*V17).*(Chi1.*V18+Chi2.*V19).^(-1);

%% Non-dimensional critical buckling loads of FGM square plates 

Ncr= Lamda0 .*a.^2./(Em*h^3);

fprintf('%.4f\n', Ncr)

