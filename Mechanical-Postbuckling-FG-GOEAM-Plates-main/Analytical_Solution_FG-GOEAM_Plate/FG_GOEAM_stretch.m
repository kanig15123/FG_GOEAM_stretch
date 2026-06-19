%% Material properties - FG-GOEAM with Quasi-3D (Thickness Stretching) Correction
T0 = 300;
DeltaT = T - T0;

EGr = 929.57e9; NuGr = 0.22; AlphaGr=-3.98e-6; RoGr=1800;
ECu = 65.79e9; NuCu = 0.387; AlphaCu=16.51e-6; RoCu=8800;

Dm = ECu*h^3/12/(1-NuCu^2);
NL = 10;

VGrs = WGr/(WGr+(RoGr/RoCu)*(1-WGr));

A11=0; A22=0; A12=0; A21=0; A66=0;
B11=0; B22=0; B12=0; B21=0; B66=0;
D11=0; D22=0; D12=0; D21=0; D66=0;
E11=0; E22=0; E12=0; E21=0; E66=0;
F11=0; F22=0; F12=0; F21=0; F66=0;
H11=0; H22=0; H12=0; H21=0; H66=0;
A44=0; A55=0; D44=0; D55=0; F44=0; F55=0;

% Quasi-3D stiffness terms
A33=0; A13=0; A23=0; NzzT=0;

% 3D Stiffness terms (un-condensed)
A11_3D=0; A12_3D=0;
B11_3D=0; B12_3D=0;
D11_3D=0; D12_3D=0;
E11_3D=0; E12_3D=0;
F11_3D=0; F12_3D=0;
H11_3D=0; H12_3D=0;

NxxT=0; NyyT=0; 
MxxT=0; MyyT=0; 
PxxT=0; PyyT=0;

for k = 1:NL
    z_bot = (NL/2-k)*h/NL; z_top = (NL/2-k+1)*h/NL;
    z_mid = (z_bot + z_top)/2;
    
    GOri_distribution;
    
    % Material modification functions
    fE = 1.11 - 1.22*VGr - 0.134*(T/T0) + 0.559*VGr*(T/T0) - 5.5*HGr*VGr + 38*HGr*VGr^2 - 20.6*HGr^2*VGr^2;
    fv = 1.01 - 1.43*VGr + 0.165*(T/T0) - 16.8*HGr*VGr - 1.1*HGr*VGr*(T/T0) + 16*HGr^2*VGr^2;
    fAlpha = 0.794 - 16.8*VGr^2 - 0.0279*(T/T0)^2 + 0.182*(T/T0)*(1+VGr);
    
    lGr=83.76e-10;
    tGr=3.4e-10;
    Xi=2*lGr/tGr;
    Eta=((EGr/ECu)-1)/((EGr/ECu)+Xi);
    VCu=1-VGr;
    
    E = (1+Xi*Eta*VGr)*ECu*fE/(1-Eta*VGr);
    Nu = (NuGr*VGr+NuCu*VCu)*fv;
    Alpha = (AlphaGr*VGr+AlphaCu*VCu)*fAlpha;
    
    % 2D Plane Stress Stiffnesses Q_ij
    constQ = 1 - Nu^2;
    Q11_val = E / constQ; 
    Q12_val = E * Nu / constQ; 
    Q66_val = E / (2*(1+Nu)); 
    Q44_val = Q66_val;
    
    % 3D Stiffnesses C_ij for Quasi-3D thickness stretching
    constC = (1 + Nu) * (1 - 2*Nu);
    C12_3D_val = E * Nu / constC; 
    C33_3D_val = E * (1 - Nu) / constC;
    
    % Define integrands (using powers of z to match HSDT formulation exactly)
    Q11a=@(z) Q11_val + 0*z; Q12a=@(z) Q12_val + 0*z; Q66a=@(z) Q66_val + 0*z;
    Q11b=@(z) Q11_val .* z; Q12b=@(z) Q12_val .* z; Q66b=@(z) Q66_val .* z;
    Q11d=@(z) Q11_val .* z.^2; Q12d=@(z) Q12_val .* z.^2; Q66d=@(z) Q66_val .* z.^2;
    Q11e=@(z) Q11_val .* z.^3; Q12e=@(z) Q12_val .* z.^3; Q66e=@(z) Q66_val .* z.^3;
    Q11f=@(z) Q11_val .* z.^4; Q12f=@(z) Q12_val .* z.^4; Q66f=@(z) Q66_val .* z.^4;
    Q11h=@(z) Q11_val .* z.^6; Q12h=@(z) Q12_val .* z.^6; Q66h=@(z) Q66_val .* z.^6;
    
    Q44a=@(z) Q44_val + 0*z; Q55a = Q44a;
    Q44d=@(z) Q44_val .* z.^2; Q55d = Q44d;
    Q44f=@(z) Q44_val .* z.^4; Q55f = Q44f;
    
    beta_val = (Q11_val + Q12_val) * Alpha * DeltaT;
    nx=@(z) beta_val + 0*z; ny=@(z) beta_val + 0*z;
    mx=@(z) beta_val .* z; my=@(z) beta_val .* z;
    px=@(z) beta_val .* z.^3; py=@(z) beta_val .* z.^3;
    
    % Thickness stretching terms
    C33_fun = @(z) C33_3D_val + 0*z;
    C13_fun = @(z) C12_3D_val + 0*z;
    Beta_z_val = (2*C12_3D_val + C33_3D_val) * Alpha * DeltaT;
    nz_fun = @(z) Beta_z_val + 0*z;
    
    % Perform integration
    nxt=integral(nx, z_bot, z_top, 'ArrayValued', true); nyt=integral(ny, z_bot, z_top, 'ArrayValued', true);
    mxt=integral(mx, z_bot, z_top, 'ArrayValued', true); myt=integral(my, z_bot, z_top, 'ArrayValued', true);
    pxt=integral(px, z_bot, z_top, 'ArrayValued', true); pyt=integral(py, z_bot, z_top, 'ArrayValued', true);
    
    a11=integral(Q11a, z_bot, z_top, 'ArrayValued', true); a12=integral(Q12a, z_bot, z_top, 'ArrayValued', true);
    a21=integral(Q12a, z_bot, z_top, 'ArrayValued', true); a22=integral(Q11a, z_bot, z_top, 'ArrayValued', true);
    a66=integral(Q66a, z_bot, z_top, 'ArrayValued', true);
    
    b11=integral(Q11b, z_bot, z_top, 'ArrayValued', true); b12=integral(Q12b, z_bot, z_top, 'ArrayValued', true);
    b21=integral(Q12b, z_bot, z_top, 'ArrayValued', true); b22=integral(Q11b, z_bot, z_top, 'ArrayValued', true);
    b66=integral(Q66b, z_bot, z_top, 'ArrayValued', true);
    
    d11=integral(Q11d, z_bot, z_top, 'ArrayValued', true); d12=integral(Q12d, z_bot, z_top, 'ArrayValued', true);
    d21=integral(Q12d, z_bot, z_top, 'ArrayValued', true); d22=integral(Q11d, z_bot, z_top, 'ArrayValued', true);
    d66=integral(Q66d, z_bot, z_top, 'ArrayValued', true);
    
    e11=integral(Q11e, z_bot, z_top, 'ArrayValued', true); e12=integral(Q12e, z_bot, z_top, 'ArrayValued', true);
    e21=integral(Q12e, z_bot, z_top, 'ArrayValued', true); e22=integral(Q11e, z_bot, z_top, 'ArrayValued', true);
    e66=integral(Q66e, z_bot, z_top, 'ArrayValued', true);
    
    f11=integral(Q11f, z_bot, z_top, 'ArrayValued', true); f12=integral(Q12f, z_bot, z_top, 'ArrayValued', true);
    f21=integral(Q12f, z_bot, z_top, 'ArrayValued', true); f22=integral(Q11f, z_bot, z_top, 'ArrayValued', true);
    f66=integral(Q66f, z_bot, z_top, 'ArrayValued', true);
    
    h11=integral(Q11h, z_bot, z_top, 'ArrayValued', true); h12=integral(Q12h, z_bot, z_top, 'ArrayValued', true);
    h21=integral(Q12h, z_bot, z_top, 'ArrayValued', true); h22=integral(Q11h, z_bot, z_top, 'ArrayValued', true);
    h66=integral(Q66h, z_bot, z_top, 'ArrayValued', true);
    
    a44=integral(Q44a, z_bot, z_top, 'ArrayValued', true); a55=a44;
    d44=integral(Q44d, z_bot, z_top, 'ArrayValued', true); d55=d44;
    f44=integral(Q44f, z_bot, z_top, 'ArrayValued', true); f55=f44;
    
    a33=integral(C33_fun, z_bot, z_top, 'ArrayValued', true);
    a13=integral(C13_fun, z_bot, z_top, 'ArrayValued', true);
    nzzt=integral(nz_fun, z_bot, z_top, 'ArrayValued', true);
    
    % Integrals for un-condensed 3D stiffnesses
    a11_3d = integral(C33_fun, z_bot, z_top, 'ArrayValued', true);
    a12_3d = integral(C13_fun, z_bot, z_top, 'ArrayValued', true);
    b11_3d = integral(@(z) C33_fun(z).*z, z_bot, z_top, 'ArrayValued', true);
    b12_3d = integral(@(z) C13_fun(z).*z, z_bot, z_top, 'ArrayValued', true);
    d11_3d = integral(@(z) C33_fun(z).*z.^2, z_bot, z_top, 'ArrayValued', true);
    d12_3d = integral(@(z) C13_fun(z).*z.^2, z_bot, z_top, 'ArrayValued', true);
    e11_3d = integral(@(z) C33_fun(z).*z.^3, z_bot, z_top, 'ArrayValued', true);
    e12_3d = integral(@(z) C13_fun(z).*z.^3, z_bot, z_top, 'ArrayValued', true);
    f11_3d = integral(@(z) C33_fun(z).*z.^4, z_bot, z_top, 'ArrayValued', true);
    f12_3d = integral(@(z) C13_fun(z).*z.^4, z_bot, z_top, 'ArrayValued', true);
    h11_3d = integral(@(z) C33_fun(z).*z.^6, z_bot, z_top, 'ArrayValued', true);
    h12_3d = integral(@(z) C13_fun(z).*z.^6, z_bot, z_top, 'ArrayValued', true);

    % Summation
    A11=A11+a11; A22=A22+a22; A12=A12+a12; A21=A21+a21; A66=A66+a66;
    B11=B11+b11; B22=B22+b22; B12=B12+b12; B21=B21+b21; B66=B66+b66;
    D11=D11+d11; D22=D22+d22; D12=D12+d12; D21=D21+d21; D66=D66+d66;
    E11=E11+e11; E22=E22+e22; E12=E12+e12; E21=E21+e21; E66=E66+e66;
    F11=F11+f11; F22=F22+f22; F12=F12+f12; F21=F21+f21; F66=F66+f66;
    H11=H11+h11; H22=H22+h22; H12=H12+h12; H21=H21+h21; H66=H66+h66;
    A44=A44+a44; A55=A55+a55; D44=D44+d44; D55=D55+d55; F44=F44+f44; F55=F55+f55;
    
    A33=A33+a33; A13=A13+a13; A23=A23+a13; NzzT=NzzT+nzzt;
    NxxT=NxxT+nxt; NyyT=NyyT+nyt; MxxT=MxxT+mxt; MyyT=MyyT+myt; PxxT=PxxT+pxt; PyyT=PyyT+pyt;
    
    A11_3D=A11_3D+a11_3d; A12_3D=A12_3D+a12_3d;
    B11_3D=B11_3D+b11_3d; B12_3D=B12_3D+b12_3d;
    D11_3D=D11_3D+d11_3d; D12_3D=D12_3D+d12_3d;
    E11_3D=E11_3D+e11_3d; E12_3D=E12_3D+e12_3d;
    F11_3D=F11_3D+f11_3d; F12_3D=F12_3D+f12_3d;
    H11_3D=H11_3D+h11_3d; H12_3D=H12_3D+h12_3d;
end