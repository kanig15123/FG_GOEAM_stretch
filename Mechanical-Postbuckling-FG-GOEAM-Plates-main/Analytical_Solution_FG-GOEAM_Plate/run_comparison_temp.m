
clear all; clc;

% Input Parameters
var1 = 15;
var2 = 1.5; 
var3 = 350;  
var4 = 85;  
var5 = 2;   
var6 = 5; 
var7 = 10;
var8 = 15; 

% Uniaxial Compression
Chi1=-1;
Chi2=0;

h=0.005;     
a=var1*h;
b=var2*a;
T=var3;
HGr=var4/100;
WGr=var5/100; 

m=1;
n=1;
c1 = 4/h^2/3;
c2 = 3*c1;

%% --- RUN WITHOUT STRETCH ---
FG_GOEAM; % Loads standard without stretch
Kun=var6; Ksn=var7; Kln=var8;
Ku=Kun*Dm/(a^4); Kl=Kln*Dm/(a^4); Ks=Ksn*Dm/(a^2);
K1=Kl*Ku/(Kl+Ku); K2=Ks*Ku/(Kl+Ku);
Parameters_Analytical_Solution;

Lambda0_ns = lamda0 * a^2 / (ECu * h^3);
Lambda1_ns = lamda1 * a^2 / (ECu * h^3);
Lambda2_ns = lamda2 * a^2 / (ECu * h^3);

%% --- RUN WITH STRETCH ---
clear Dm ECu h a b T HGr WGr m n c1 c2 K1 K2 Ku Kl Ks lamda0 lamda1 lamda2;
h=0.005;     
a=var1*h;
b=var2*a;
T=var3;
HGr=var4/100;
WGr=var5/100; 
m=1;
n=1;
c1 = 4/h^2/3;
c2 = 3*c1;

FG_GOEAM_stretch; % Loads stretch
Kun=var6; Ksn=var7; Kln=var8;
Ku=Kun*Dm/(a^4); Kl=Kln*Dm/(a^4); Ks=Ksn*Dm/(a^2);
K1=Kl*Ku/(Kl+Ku); K2=Ks*Ku/(Kl+Ku);
Parameters_stretch;

Lambda0_s = lamda0 * a^2 / (ECu * h^3);
Lambda1_s = lamda1 * a^2 / (ECu * h^3);
Lambda2_s = lamda2 * a^2 / (ECu * h^3);

%% --- PRINT RESULTS ---
fprintf('RESULT_START\n');
fprintf('No-Stretch: Lambda0=%.6f, Lambda1=%.6f, Lambda2=%.6f\n', Lambda0_ns, Lambda1_ns, Lambda2_ns);
fprintf('Stretch:    Lambda0=%.6f, Lambda1=%.6f, Lambda2=%.6f\n', Lambda0_s, Lambda1_s, Lambda2_s);
fprintf('RESULT_END\n');

%% --- PLOT AND SAVE ---
num_points = 100; 
Wn = linspace(0, 1, num_points);
Lambda_n_ns = Lambda0_ns + Lambda1_ns .* Wn + Lambda2_ns .* (Wn.^2);
Lambda_n_s = Lambda0_s + Lambda1_s .* Wn + Lambda2_s .* (Wn.^2);

fig = figure('Visible', 'off');
plot(Wn, Lambda_n_ns, 'k--', 'LineWidth', 2.5); hold on;
plot(Wn, Lambda_n_s, 'r-', 'LineWidth', 2.5);
grid on; box on;
set(gca, 'LineWidth', 1.5, 'FontSize', 14);
xlabel('\overline{W}');
ylabel('\overline{\lambda}');
legend('Standard HSDT (No Stretch)', 'Quasi-3D (With Stretch)', 'Location', 'Best');
title('Postbuckling Load-Deflection Curve Comparison');
saveas(fig, 'postbuckling_comparison.png');
