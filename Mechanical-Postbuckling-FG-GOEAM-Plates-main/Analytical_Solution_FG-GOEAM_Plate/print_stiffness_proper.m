
clear all; clc;

var1 = 15; var2 = 1.5; var3 = 350; var4 = 85; var5 = 2; var6 = 5; var7 = 10; var8 = 15; 
Chi1=-1; Chi2=0;
h=0.005; a=var1*h; b=var2*a; T=var3; HGr=var4/100; WGr=var5/100; 
m=1; n=1; c1 = 4/h^2/3; c2 = 3*c1;

%% --- RUN WITHOUT STRETCH ---
FG_GOEAM;
Kun=var6; Ksn=var7; Kln=var8;
Ku=Kun*Dm/(a^4); Kl=Kln*Dm/(a^4); Ks=Ksn*Dm/(a^2);
K1=Kl*Ku/(Kl+Ku); K2=Ks*Ku/(Kl+Ku);

fprintf('=== NO STRETCH ===\n');
fprintf('A11 = %.4e, A12 = %.4e, A22 = %.4e, A66 = %.4e\n', A11, A12, A22, A66);
fprintf('B11 = %.4e, B12 = %.4e, B22 = %.4e, B66 = %.4e\n', B11, B12, B22, B66);
fprintf('D11 = %.4e, D12 = %.4e, D22 = %.4e, D66 = %.4e\n', D11, D12, D22, D66);
Parameters_Analytical_Solution;
fprintf('Denom = %.4e\n', (A12*A21 - A11*A22));
fprintf('C11 = %.4e, C12 = %.4e, C13 = %.4e, C14 = %.4e\n', C11, C12, C13, C14);
fprintf('J11 = %.4e, J12 = %.4e, J13 = %.4e\n', J11, J12, J13);
fprintf('R11 = %.4e, R21 = %.4e, R31 = %.4e, R32 = %.4e, R33 = %.4e\n', R11, R21, R31, R32, R33);
fprintf('lamda0 = %.4e, lamda1 = %.4e, lamda2 = %.4e\n', lamda0, lamda1, lamda2);

%% --- RUN WITH CORRECTED PROPER STRETCH ---
clear Dm ECu h a b T HGr WGr m n c1 c2 K1 K2 Ku Kl Ks lamda0 lamda1 lamda2;
h=0.005; a=var1*h; b=var2*a; T=var3; HGr=var4/100; WGr=var5/100; 
m=1; n=1; c1 = 4/h^2/3; c2 = 3*c1;

T0 = 300; DeltaT = T - T0; 
EGr = 929.57e9; NuGr = 0.22; AlphaGr=-3.98e-6; RoGr=1800;
ECu = 65.79e9; NuCu = 0.387; AlphaCu=16.51e-6; RoCu=8800;
Dm = ECu * h^3 / 12 / (1 - NuCu^2);
NL = 10; dz = h / NL;
VGrs = WGr / (WGr + (RoGr/RoCu)*(1 - WGr));

A11=0; A22=0; A12=0; A21=0; A66=0;
B11=0; B22=0; B12=0; B21=0; B66=0;
D11=0; D22=0; D12=0; D21=0; D66=0;
E11=0; E22=0; E12=0; E21=0; E66=0;
F11=0; F22=0; F12=0; F21=0; F66=0;
H11=0; H22=0; H12=0; H21=0; H66=0;
A44=0; A55=0; D44=0; D55=0; F44=0; F55=0;
A33=0; A13=0; A23=0; NzzT=0;
NxxT=0; NyyT=0; MxxT=0; MyyT=0; PxxT=0; PyyT=0;

for k = 1:NL
    z_bot = (NL/2-k)*h/NL; z_top = (NL/2-k+1)*h/NL;
    GOri_distribution; 
    TR = T / T0;
    fE = 1.11 - 1.22*VGr - 0.134*TR + 0.559*VGr*TR - 5.5*HGr*VGr + 38*HGr*VGr^2 - 20.6*HGr^2*VGr^2;
    fv = 1.01 - 1.43*VGr + 0.165*TR - 16.8*HGr*VGr - 1.1*HGr*VGr*TR + 16*HGr^2*VGr^2;
    fAlpha = 0.794 - 16.8*VGr^2 - 0.0279*TR^2 + 0.182*TR*(1+VGr);
    
    lGr=83.76e-10; tGr=3.4e-10;
    Xi=2*lGr/tGr;
    Eta=((EGr/ECu)-1)/((EGr/ECu)+Xi);
    VCu=1-VGr;
    
    E = (1+Xi*Eta*VGr)*ECu*fE/(1-Eta*VGr);
    Nu = (NuGr*VGr+NuCu*VCu)*fv;
    Alpha = (AlphaGr*VGr+AlphaCu*VCu)*fAlpha;

    f_fun = @(z) z - c1*z.^3; g_fun = @(z) 1 - 3*c1*z.^2;
    
    constQ = 1 - Nu^2;
    Q11_val = E / constQ; Q12_val = E * Nu / constQ; Q66_val = E / (2*(1+Nu)); Q44_val = Q66_val; 
    
    constC = (1 + Nu) * (1 - 2*Nu);
    C11_3D_val = E * (1 - Nu) / constC; C12_3D_val = E * Nu / constC; C33_3D_val = C11_3D_val; 
    
    Q11a = @(z) Q11_val + 0*z; Q12a = @(z) Q12_val + 0*z; Q66a = @(z) Q66_val + 0*z;
    Q11b = @(z) Q11_val .* z; Q12b = @(z) Q12_val .* z; Q66b = @(z) Q66_val .* z;
    Q11d = @(z) Q11_val .* z.^2; Q12d = @(z) Q12_val .* z.^2; Q66d = @(z) Q66_val .* z.^2;
    Q11e = @(z) Q11_val .* f_fun(z); Q12e = @(z) Q12_val .* f_fun(z); Q66e = @(z) Q66_val .* f_fun(z);
    Q11f = @(z) Q11_val .* z .* f_fun(z); Q12f = @(z) Q12_val .* z .* f_fun(z); Q66f = @(z) Q66_val .* z .* f_fun(z);
    Q11h = @(z) Q11_val .* f_fun(z).^2; Q12h = @(z) Q12_val .* f_fun(z).^2; Q66h = @(z) Q66_val .* f_fun(z).^2;
    Q44a = @(z) Q44_val * (g_fun(z).^2); Q55a = Q44a; 
    Q44d = @(z) Q44_val * (g_fun(z).^2) .* z.^2; Q55d = Q44d;
    Q44f = @(z) Q44_val * (g_fun(z).^2) .* f_fun(z); Q55f = Q44f;
    
    beta_val = (Q11_val + Q12_val) * Alpha * DeltaT;
    nx = @(z) beta_val + 0*z; ny = nx; mx = @(z) beta_val .* z; my = mx; px = @(z) beta_val .* f_fun(z); py = px;
    
    C33_fun = @(z) C33_3D_val + 0*z; C13_fun = @(z) C12_3D_val + 0*z;
    Beta_z_val = (2*C12_3D_val + C33_3D_val) * Alpha * DeltaT;
    nz_fun = @(z) Beta_z_val + 0*z;

    nxt=integral(nx, z_bot, z_top, 'ArrayValued', true); nyt=integral(ny, z_bot, z_top, 'ArrayValued', true);
    mxt=integral(mx, z_bot, z_top, 'ArrayValued', true); myt=integral(my, z_bot, z_top, 'ArrayValued', true);
    pxt=integral(px, z_bot, z_top, 'ArrayValued', true); pyt=integral(py, z_bot, z_top, 'ArrayValued', true);
    
    a11=integral(Q11a, z_bot, z_top, 'ArrayValued', true); 
    a12=integral(Q12a, z_bot, z_top, 'ArrayValued', true);
    a21=integral(Q12a, z_bot, z_top, 'ArrayValued', true); 
    a22=integral(Q11a, z_bot, z_top, 'ArrayValued', true); 
    a66=integral(Q66a, z_bot, z_top, 'ArrayValued', true);
    
    b11=integral(Q11b, z_bot, z_top, 'ArrayValued', true); 
    b12=integral(Q12b, z_bot, z_top, 'ArrayValued', true);
    b21=integral(Q12b, z_bot, z_top, 'ArrayValued', true); 
    b22=integral(Q11b, z_bot, z_top, 'ArrayValued', true); 
    b66=integral(Q66b, z_bot, z_top, 'ArrayValued', true);
    
    d11=integral(Q11d, z_bot, z_top, 'ArrayValued', true); 
    d12=integral(Q12d, z_bot, z_top, 'ArrayValued', true);
    d21=integral(Q12d, z_bot, z_top, 'ArrayValued', true); 
    d22=integral(Q11d, z_bot, z_top, 'ArrayValued', true); 
    d66=integral(Q66d, z_bot, z_top, 'ArrayValued', true);
    
    e11=integral(Q11e, z_bot, z_top, 'ArrayValued', true); 
    e12=integral(Q12e, z_bot, z_top, 'ArrayValued', true);
    e21=integral(Q12e, z_bot, z_top, 'ArrayValued', true); 
    e22=integral(Q11e, z_bot, z_top, 'ArrayValued', true); 
    e66=integral(Q66e, z_bot, z_top, 'ArrayValued', true);
    
    f11=integral(Q11f, z_bot, z_top, 'ArrayValued', true); 
    f12=integral(Q12f, z_bot, z_top, 'ArrayValued', true);
    f21=integral(Q12f, z_bot, z_top, 'ArrayValued', true); 
    f22=integral(Q11f, z_bot, z_top, 'ArrayValued', true); 
    f66=integral(Q66f, z_bot, z_top, 'ArrayValued', true);
    
    h11=integral(Q11h, z_bot, z_top, 'ArrayValued', true); 
    h12=integral(Q12h, z_bot, z_top, 'ArrayValued', true);
    h21=integral(Q12h, z_bot, z_top, 'ArrayValued', true); 
    h22=integral(Q11h, z_bot, z_top, 'ArrayValued', true); 
    h66=integral(Q66h, z_bot, z_top, 'ArrayValued', true);
    
    a44=integral(Q44a, z_bot, z_top, 'ArrayValued', true); a55=a44;
    d44=integral(Q44d, z_bot, z_top, 'ArrayValued', true); d55=d44;
    f44=integral(Q44f, z_bot, z_top, 'ArrayValued', true); f55=f44;
    
    a33 = integral(C33_fun, z_bot, z_top, 'ArrayValued', true);
    a13 = integral(C13_fun, z_bot, z_top, 'ArrayValued', true);
    a23 = a13;
    nzzt = integral(nz_fun, z_bot, z_top, 'ArrayValued', true);

    A11=A11+a11; A22=A22+a22; A12=A12+a12; A21=A21+a21; A66=A66+a66;
    B11=B11+b11; B22=B22+b22; B12=B12+b12; B21=B21+b21; B66=B66+b66;
    D11=D11+d11; D22=D22+d22; D12=D12+d12; D21=D21+d21; D66=D66+d66;
    E11=E11+e11; E22=E22+e22; E12=E12+e12; E21=E21+e21; E66=E66+e66;
    F11=F11+f11; F22=F22+f22; F12=F12+f12; F21=F21+f21; F66=F66+f66;
    H11=H11+h11; H22=H22+h22; H12=H12+h12; H21=H21+h21; H66=H66+h66;
    A44=A44+a44; A55=A55+a55; D44=D44+d44; D55=D55+d55; F44=F44+f44; F55=F55+f55;
    A33=A33+a33; A13=A13+a13; A23=A23+a23; NzzT=NzzT+nzzt;
    NxxT=NxxT+nxt; NyyT=NyyT+nyt; MxxT=MxxT+mxt; MyyT=MyyT+myt; PxxT=PxxT+pxt; PyyT=PyyT+pyt;
end

fprintf('=== WITH STRETCH (BEFORE CONDENSATION) ===\n');
fprintf('A11 = %.4e, A12 = %.4e, A22 = %.4e, A66 = %.4e\n', A11, A12, A22, A66);
fprintf('A13 = %.4e, A23 = %.4e, A33 = %.4e\n', A13, A23, A33);

Kun=var6; Ksn=var7; Kln=var8;
Ku=Kun*Dm/(a^4); Kl=Kln*Dm/(a^4); Ks=Ksn*Dm/(a^2);
K1=Kl*Ku/(Kl+Ku); K2=Ks*Ku/(Kl+Ku);

Parameters_stretch;

fprintf('=== WITH STRETCH (AFTER CONDENSATION) ===\n');
fprintf('A11 = %.4e, A12 = %.4e, A22 = %.4e, A66 = %.4e\n', A11, A12, A22, A66);
fprintf('Denom = %.4e\n', (A12*A21 - A11*A22));
fprintf('C11 = %.4e, C12 = %.4e, C13 = %.4e, C14 = %.4e\n', C11, C12, C13, C14);
fprintf('J11 = %.4e, J12 = %.4e, J13 = %.4e\n', J11, J12, J13);
fprintf('R11 = %.4e, R21 = %.4e, R31 = %.4e, R32 = %.4e, R33 = %.4e\n', R11, R21, R31, R32, R33);
fprintf('lamda0 = %.4e, lamda1 = %.4e, lamda2 = %.4e\n', lamda0, lamda1, lamda2);

