%% ========================================================================
%  Project:       Mechanical Postbuckling of FG-GOEAM Plates
%  Description:   Comparative Analysis of Standard HSDT vs. Quasi-3D HSDT
%                 (incorporating the Thickness Stretching effect)
% =========================================================================

clear all; clc; close all;

%% ------------------------- Input Parameters ----------------------------
var1 = 15;   % Aspect ratio (a/h)
var2 = 1.5;  % Length-to-width ratio (b/a)
var3 = 350;  % Temperature T (K)
var4 = 85;   % Degree of GOri folding HGr (%)
var5 = 2;    % GOri weight fraction WGr (%)
var6 = 5;    % Winkler spring stiffness (Kun)
var7 = 10;   % Pasternak shear layer stiffness (Ksn)
var8 = 15;   % Lower spring stiffness (Ksn)

h = 0.005;   % Thickness (m)
a = var1*h;  % Length (m)
b = var2*a;  % Width (m)
T = var3;    
HGr = var4/100;
WGr = var5/100;

% Mode numbers
m = 1;
n = 1;

% Constants for HSDT formulation
c1 = 4/h^2/3;
c2 = 3*c1;

%% ------------------------- Loading Case --------------------------------
% Uniaxial Compression
Chi1 = -1;
Chi2 = 0;

% Biaxial Compression (uncomment to test)
% Chi1 = -1;
% Chi2 = -1;

%% ========================================================================
%  1. RUN WITHOUT THICKNESS STRETCHING (BASELINE HSDT)
% =========================================================================
FG_GOEAM; % Load standard HSDT properties

Kun = var6; Ksn = var7; Kln = var8;
Ku = Kun*Dm/(a^4); Kl = Kln*Dm/(a^4); Ks = Ksn*Dm/(a^2);
K1 = Kl*Ku/(Kl+Ku); K2 = Ks*Ku/(Kl+Ku);

Parameters_Analytical_Solution; % Solve baseline parameters

Lambda0_ns = lamda0 * a^2 / (ECu * h^3);
Lambda1_ns = lamda1 * a^2 / (ECu * h^3);
Lambda2_ns = lamda2 * a^2 / (ECu * h^3);

%% ========================================================================
%  2. RUN WITH THICKNESS STRETCHING (QUASI-3D HSDT)
% ========================================================================
% Clear intermediate variables to avoid pollution
clear Dm ECu h a b T HGr WGr m n c1 c2 K1 K2 Ku Kl Ks lamda0 lamda1 lamda2;

h = 0.005;
a = var1*h;
b = var2*a;
T = var3;
HGr = var4/100;
WGr = var5/100;
m = 1; n = 1;
c1 = 4/h^2/3; c2 = 3*c1;

FG_GOEAM_stretch; % Load properties with thickness stretching terms

Kun = var6; Ksn = var7; Kln = var8;
Ku = Kun*Dm/(a^4); Kl = Kln*Dm/(a^4); Ks = Ksn*Dm/(a^2);
K1 = Kl*Ku/(Kl+Ku); K2 = Ks*Ku/(Kl+Ku);

Parameters_stretch; % Solve with static condensation applied

Lambda0_s = lamda0 * a^2 / (ECu * h^3);
Lambda1_s = lamda1 * a^2 / (ECu * h^3);
Lambda2_s = lamda2 * a^2 / (ECu * h^3);

%% ========================================================================
%  3. DISPLAY NUMERICAL RESULTS
% ========================================================================
fprintf('\n========================================================================\n');
fprintf('                POSTBUCKLING COEFFICIENTS COMPARISON\n');
fprintf('========================================================================\n');
fprintf('Coefficient    No-Stretch (HSDT)    Stretch (Quasi-3D)    %% Change\n');
fprintf('------------------------------------------------------------------------\n');
fprintf('Lambda0        %17.6f    %18.6f    %7.4f%%\n', ...
    Lambda0_ns, Lambda0_s, ((Lambda0_ns - Lambda0_s)/Lambda0_ns)*100);
fprintf('Lambda1        %17.6f    %18.6f    %7.4f%%\n', ...
    Lambda1_ns, Lambda1_s, 0.0); % Avoid division by zero for lambda1 ~ 0
fprintf('Lambda2        %17.6f    %18.6f    %7.4f%%\n', ...
    Lambda2_ns, Lambda2_s, ((Lambda2_ns - Lambda2_s)/Lambda2_ns)*100);
fprintf('========================================================================\n\n');

%% ========================================================================
%  4. PLOT LOAD-DEFLECTION CURVES
% ========================================================================
num_points = 100;
Wn = linspace(0, 1, num_points); % Normalized deflection W_bar

Lambda_n_ns = Lambda0_ns + Lambda1_ns .* Wn + Lambda2_ns .* (Wn.^2);
Lambda_n_s  = Lambda0_s  + Lambda1_s  .* Wn + Lambda2_s  .* (Wn.^2);

figure;
plot(Wn, Lambda_n_ns, 'k--', 'LineWidth', 3.0); hold on;
plot(Wn, Lambda_n_s, 'r-', 'LineWidth', 3.0);

grid on; box on;
set(gca, 'LineWidth', 2, 'FontSize', 22, 'TickLabelInterpreter', 'latex');
pbaspect([4 3 1]);

xlabel('$\overline{W}$', 'FontSize', 26, 'Interpreter', 'latex');
ylabel('$\overline{\lambda}$', 'FontSize', 26, 'Interpreter', 'latex');
title('\textbf{Postbuckling Comparison: Standard HSDT vs. Quasi-3D}', ...
      'FontSize', 24, 'FontWeight', 'bold', 'Interpreter', 'latex');
legend({'Standard HSDT (No Stretch)', 'Quasi-3D HSDT (With Stretch)'}, ...
       'FontSize', 18, 'Interpreter', 'latex', 'Location', 'NorthWest');

% Save comparison figure
saveas(gcf, 'postbuckling_comparison_corrected.png');
fprintf('Figure saved as "postbuckling_comparison_corrected.png"\n');
