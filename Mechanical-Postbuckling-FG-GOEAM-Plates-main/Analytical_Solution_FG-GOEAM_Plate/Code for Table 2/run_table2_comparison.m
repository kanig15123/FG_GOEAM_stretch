
clear all; clc;

%% --- RUN WITHOUT STRETCH ---
h=0.002; a=10*h; b=a; m=1; n=1; c1 = 4/h^2/3; c2 = 3*c1; p = 1;
FGM; 
Dm=Em*h^3/12/(1-Num^2);
Kwn=100; Gpn=10;
Kw=Kwn*Dm/a^4; Gp=Gpn*Dm/a^2;
K1=Kw; K2=Gp;
Chi1=-1; Chi2=-1;

Parameters_Analytical_Solution;
Lamda0=(-1).*(S2+V11+P11.*V14+P21.*V15).*(Chi1.*V18+Chi2.*V19).^(-1);
Ncr_ns = Lamda0 * a^2 / (Em*h^3);

%% --- RUN WITH STRETCH (STATIC CONDENSATION) ---
clear Dm Kw Gp K1 K2 lamda0 lamda1 lamda2;
h=0.002; a=10*h; b=a; m=1; n=1; c1 = 4/h^2/3; c2 = 3*c1; p = 1;

% Load properties but add thickness stretching stiffnesses
T=300; T0=300; DeltaT=T-T0;
Ec = 380E9; Nuc = 0.3; Alphac = 7.4e-6;
Em = 70E9; Num = 0.3; Alpham = 23e-6;

E = @(z) Em + (Ec - Em) .* (h.^(-1) .* ((1./2) * h + z)).^p;
Nu = @(z) Num + (Nuc - Num) .* (h.^(-1) .* ((1./2) .* h + z)).^p;
Alpha = @(z) Alpham + (Alphac - Alpham) .* (h.^(-1) .* ((1./2) .* h + z)).^p;

Q11a=@(z) (E(z)./(1-Nu(z).^2)); Q12a=@(z) (E(z).*Nu(z)./(1-Nu(z).^2)); Q66a=@(z) (E(z)./2./(1+Nu(z)));
Q11b=@(z) (E(z)./(1-Nu(z).^2)).*z; Q12b=@(z) (E(z).*Nu(z)./(1-Nu(z).^2)).*z; Q66b=@(z) (E(z)./2./(1+Nu(z))).*z;
Q11d=@(z) (E(z)./(1-Nu(z).^2)).*(z.^2); Q12d=@(z) (E(z).*Nu(z)./(1-Nu(z).^2)).*(z.^2); Q66d=@(z) (E(z)./2./(1+Nu(z))).*(z.^2);
Q11e=@(z) (E(z)./(1-Nu(z).^2)).*(z.^3); Q12e=@(z) (E(z).*Nu(z)./(1-Nu(z).^2)).*(z.^3); Q66e=@(z) (E(z)./2./(1+Nu(z))).*(z.^3);
Q11f=@(z) (E(z)./(1-Nu(z).^2)).*(z.^4); Q12f=@(z) (E(z).*Nu(z)./(1-Nu(z).^2)).*(z.^4); Q66f=@(z) (E(z)./2./(1+Nu(z))).*(z.^4);
Q11h=@(z) (E(z)./(1-Nu(z).^2)).*(z.^6); Q12h=@(z) (E(z).*Nu(z)./(1-Nu(z).^2)).*(z.^6); Q66h=@(z) (E(z)./2./(1+Nu(z))).*(z.^6);
Q44a=@(z) (E(z)./2./(1+Nu(z))); Q55a=Q44a;
Q44d=@(z) (E(z)./2./(1+Nu(z))).*(z.^2); Q55d=Q44d;
Q44f=@(z) (E(z)./2./(1+Nu(z))).*(z.^4); Q55f=Q44f;

constC = @(z) (1 + Nu(z)) .* (1 - 2*Nu(z));
C33_fun = @(z) E(z) .* (1 - Nu(z)) ./ constC(z);
C13_fun = @(z) E(z) .* Nu(z) ./ constC(z);

nx=@(z) (Q11a(z)+Q12a(z)).*Alpha(z).*DeltaT; ny=@(z) (Q11a(z)+Q12a(z)).*Alpha(z).*DeltaT;
mx=@(z) (Q11a(z)+Q12a(z)).*Alpha(z).*DeltaT.*z; my=@(z) (Q11a(z)+Q12a(z)).*Alpha(z).*DeltaT.*z;
px=@(z) (Q11a(z)+Q12a(z)).*Alpha(z).*DeltaT.*(z.^3); py=@(z) (Q11a(z)+Q12a(z)).*Alpha(z).*DeltaT.*(z.^3);
nz_fun = @(z) (2*C13_fun(z) + C33_fun(z)) .* Alpha(z) .* DeltaT;

NxxT=integral(nx, -h/2, h/2); NyyT=integral(ny, -h/2, h/2); MxxT=integral(mx, -h/2, h/2); MyyT=integral(my, -h/2, h/2); PxxT=integral(px, -h/2, h/2); PyyT=integral(py, -h/2, h/2);
A11=integral(Q11a, -h/2, h/2); A12=integral(Q12a, -h/2, h/2); A21=A12; A22=A11; A66=integral(Q66a, -h/2, h/2);
B11=integral(Q11b, -h/2, h/2); B12=integral(Q12b, -h/2, h/2); B21=B12; B22=B11; B66=integral(Q66b, -h/2, h/2);
D11=integral(Q11d, -h/2, h/2); D12=integral(Q12d, -h/2, h/2); D21=D12; D22=D11; D66=integral(Q66d, -h/2, h/2);
E11=integral(Q11e, -h/2, h/2); E12=integral(Q12e, -h/2, h/2); E21=E12; E22=E11; E66=integral(Q66e, -h/2, h/2);
F11=integral(Q11f, -h/2, h/2); F12=integral(Q12f, -h/2, h/2); F21=F12; F22=F11; F66=integral(Q66f, -h/2, h/2);
H11=integral(Q11h, -h/2, h/2); H12=integral(Q12h, -h/2, h/2); H21=H12; H22=H11; H66=integral(Q66h, -h/2, h/2);
A44=integral(Q44a, -h/2, h/2); A55=A44; D44=integral(Q44d, -h/2, h/2); D55=D44; F44=integral(Q44f, -h/2, h/2); F55=F44;

A33 = integral(C33_fun, -h/2, h/2);
A13 = integral(C13_fun, -h/2, h/2);
A23 = A13;
NzzT = integral(nz_fun, -h/2, h/2);

Dm=Em*h^3/12/(1-Num^2);
Kwn=100; Gpn=10; Kw=Kwn*Dm/a^4; Gp=Gpn*Dm/a^2; K1=Kw; K2=Gp;
Chi1=-1; Chi2=-1;

% Call parameters_stretch to apply condensation (from parent folder)
addpath('..');
Parameters_stretch;
Ncr_s = lamda0 * a^2 / (Em*h^3);

fprintf('TABLE2_START\n');
fprintf('Table 2 - No-Stretch: %.6f\n', Ncr_ns);
fprintf('Table 2 - Stretch:    %.6f\n', Ncr_s);
fprintf('Table 2 - Diff:       %.6f\n', Ncr_ns - Ncr_s);
fprintf('TABLE2_END\n');


