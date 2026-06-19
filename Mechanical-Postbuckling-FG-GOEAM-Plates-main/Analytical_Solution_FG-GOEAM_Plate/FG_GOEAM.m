%% Material properties- FG-GOEAM
T0=300;
DeltaT=T-T0;

EGr = 929.57e9; NuGr = 0.22; AlphaGr=-3.98e-6; RoGr=1800;
ECu = 65.79e9; NuCu = 0.387; AlphaCu=16.51e-6; RoCu=8800;

Dm=ECu*h^3/12/(1-NuCu^2);

NL=10;

VGrs=WGr/(WGr+(RoGr/RoCu)*(1-WGr));

A11=0; A22=0; A12=0; A21=0; A66=0;
B11=0; B22=0; B12=0; B21=0; B66=0;
D11=0; D22=0; D12=0; D21=0; D66=0;
E11=0; E22=0; E12=0; E21=0; E66=0;
F11=0; F22=0; F12=0; F21=0; F66=0;
H11=0; H22=0; H12=0; H21=0; H66=0;
A44=0; A55=0; D44=0; D55=0; F44=0; 
F55=0;

NxxT=0; NyyT=0; 
MxxT=0; MyyT=0; 
PxxT=0; PyyT=0;

for k=1:1:NL
    
GOri_distribution

% material functions
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

% Qij matrix
Q11a=@(z) (E/(1-Nu^2)); Q22a=@(z) (E/(1-Nu^2)); Q12a=@(z) (E*Nu/(1-Nu^2));
Q21a=@(z) (E*Nu/(1-Nu^2)); Q66a=@(z) (E/2/(1+Nu));

Q11b=@(z) (E/(1-Nu^2)).*z; Q22b=@(z) (E/(1-Nu^2)).*z; Q12b=@(z) (E*Nu/(1-Nu^2)).*z;
Q21b=@(z) (E*Nu/(1-Nu^2)).*z; Q66b=@(z) (E/2/(1+Nu)).*z;

Q11d=@(z) (E/(1-Nu^2)).*(z.^2); Q22d=@(z) (E/(1-Nu^2)).*(z.^2); Q12d=@(z) (E*Nu/(1-Nu^2)).*(z.^2);
Q21d=@(z) (E*Nu/(1-Nu^2)).*(z.^2); Q66d=@(z) (E/2/(1+Nu)).*(z.^2);

Q11e=@(z) (E/(1-Nu^2)).*(z.^3); Q22e=@(z) (E/(1-Nu^2)).*(z.^3); Q12e=@(z) (E*Nu/(1-Nu^2)).*(z.^3);
Q21e=@(z) (E*Nu/(1-Nu^2)).*(z.^3); Q66e=@(z) (E/2/(1+Nu)).*(z.^3);

Q11f=@(z) (E/(1-Nu^2)).*(z.^4); Q22f=@(z) (E/(1-Nu^2)).*(z.^4); Q12f=@(z) (E*Nu/(1-Nu^2)).*(z.^4);
Q21f=@(z) (E*Nu/(1-Nu^2)).*(z.^4); Q66f=@(z) (E/2/(1+Nu)).*(z.^4);

Q11h=@(z) (E/(1-Nu^2)).*(z.^6); Q22h=@(z) (E/(1-Nu^2)).*(z.^6); Q12h=@(z) (E*Nu/(1-Nu^2)).*(z.^6);
Q21h=@(z) (E*Nu/(1-Nu^2)).*(z.^6); Q66h=@(z) (E/2/(1+Nu)).*(z.^6);

Q44a=@(z) (E/2/(1+Nu)); Q55a=@(z) (E/2/(1+Nu));
Q44d=@(z) (E/2/(1+Nu)).*(z.^2); Q55d=@(z) (E/2/(1+Nu)).*(z.^2);
Q44f=@(z) (E/2/(1+Nu)).*(z.^4); Q55f=@(z) (E/2/(1+Nu)).*(z.^4);


nx=@(z) (Q11a(z)+Q12a(z)).*Alpha.*DeltaT; ny=@(z) (Q21a(z)+Q22a(z)).*Alpha.*DeltaT;
mx=@(z) (Q11a(z)+Q12a(z)).*Alpha.*DeltaT.*z; my=@(z) (Q21a(z)+Q22a(z)).*Alpha.*DeltaT.*z;
px=@(z) (Q11a(z)+Q12a(z)).*Alpha.*DeltaT.*(z.^3); py=@(z) (Q21a(z)+Q22a(z)).*Alpha.*DeltaT.*(z.^3);


%%%%%%%%%%%%%%%%%  plate stiffness coefficients %%%%%%%%%%%%%%%%%%%%%%%%%%%
nxt=integral(nx,(NL/2-k)*h/NL,(NL/2-k+1)*h/NL,'ArrayValued',true); nyt=integral(ny,(NL/2-k)*h/NL,(NL/2-k+1)*h/NL,'ArrayValued',true);
mxt=integral(mx,(NL/2-k)*h/NL,(NL/2-k+1)*h/NL,'ArrayValued',true); myt=integral(my,(NL/2-k)*h/NL,(NL/2-k+1)*h/NL,'ArrayValued',true);
pxt=integral(px,(NL/2-k)*h/NL,(NL/2-k+1)*h/NL,'ArrayValued',true); pyt=integral(py,(NL/2-k)*h/NL,(NL/2-k+1)*h/NL,'ArrayValued',true);

a11=integral(Q11a,(NL/2-k)*h/NL,(NL/2-k+1)*h/NL,'ArrayValued',true); a12=integral(Q12a,(NL/2-k)*h/NL,(NL/2-k+1)*h/NL,'ArrayValued',true);
a21=integral(Q21a,(NL/2-k)*h/NL,(NL/2-k+1)*h/NL,'ArrayValued',true); a22=integral(Q22a,(NL/2-k)*h/NL,(NL/2-k+1)*h/NL,'ArrayValued',true);
a66=integral(Q66a,(NL/2-k)*h/NL,(NL/2-k+1)*h/NL,'ArrayValued',true);

b11=integral(Q11b,(NL/2-k)*h/NL,(NL/2-k+1)*h/NL,'ArrayValued',true); b12=integral(Q12b,(NL/2-k)*h/NL,(NL/2-k+1)*h/NL,'ArrayValued',true);
b21=integral(Q21b,(NL/2-k)*h/NL,(NL/2-k+1)*h/NL,'ArrayValued',true); b22=integral(Q22b,(NL/2-k)*h/NL,(NL/2-k+1)*h/NL,'ArrayValued',true);
b66=integral(Q66b,(NL/2-k)*h/NL,(NL/2-k+1)*h/NL,'ArrayValued',true);

d11=integral(Q11d,(NL/2-k)*h/NL,(NL/2-k+1)*h/NL,'ArrayValued',true); d12=integral(Q12d,(NL/2-k)*h/NL,(NL/2-k+1)*h/NL,'ArrayValued',true);
d21=integral(Q21d,(NL/2-k)*h/NL,(NL/2-k+1)*h/NL,'ArrayValued',true); d22=integral(Q22d,(NL/2-k)*h/NL,(NL/2-k+1)*h/NL,'ArrayValued',true);
d66=integral(Q66d,(NL/2-k)*h/NL,(NL/2-k+1)*h/NL,'ArrayValued',true);

e11=integral(Q11e,(NL/2-k)*h/NL,(NL/2-k+1)*h/NL,'ArrayValued',true); e12=integral(Q12e,(NL/2-k)*h/NL,(NL/2-k+1)*h/NL,'ArrayValued',true);
e21=integral(Q21e,(NL/2-k)*h/NL,(NL/2-k+1)*h/NL,'ArrayValued',true); e22=integral(Q22e,(NL/2-k)*h/NL,(NL/2-k+1)*h/NL,'ArrayValued',true);
e66=integral(Q66e,(NL/2-k)*h/NL,(NL/2-k+1)*h/NL,'ArrayValued',true);

f11=integral(Q11f,(NL/2-k)*h/NL,(NL/2-k+1)*h/NL,'ArrayValued',true); f12=integral(Q12f,(NL/2-k)*h/NL,(NL/2-k+1)*h/NL,'ArrayValued',true);
f21=integral(Q21f,(NL/2-k)*h/NL,(NL/2-k+1)*h/NL,'ArrayValued',true); f22=integral(Q22f,(NL/2-k)*h/NL,(NL/2-k+1)*h/NL,'ArrayValued',true);
f66=integral(Q66f,(NL/2-k)*h/NL,(NL/2-k+1)*h/NL,'ArrayValued',true);

h11=integral(Q11h,(NL/2-k)*h/NL,(NL/2-k+1)*h/NL,'ArrayValued',true); h12=integral(Q12h,(NL/2-k)*h/NL,(NL/2-k+1)*h/NL,'ArrayValued',true);
h21=integral(Q21h,(NL/2-k)*h/NL,(NL/2-k+1)*h/NL,'ArrayValued',true); h22=integral(Q22h,(NL/2-k)*h/NL,(NL/2-k+1)*h/NL,'ArrayValued',true);
h66=integral(Q66h,(NL/2-k)*h/NL,(NL/2-k+1)*h/NL,'ArrayValued',true);

a44=integral(Q44a,(NL/2-k)*h/NL,(NL/2-k+1)*h/NL,'ArrayValued',true); a55=integral(Q55a,(NL/2-k)*h/NL,(NL/2-k+1)*h/NL,'ArrayValued',true);
d44=integral(Q44d,(NL/2-k)*h/NL,(NL/2-k+1)*h/NL,'ArrayValued',true); d55=integral(Q55d,(NL/2-k)*h/NL,(NL/2-k+1)*h/NL,'ArrayValued',true);
f44=integral(Q44f,(NL/2-k)*h/NL,(NL/2-k+1)*h/NL,'ArrayValued',true); f55=integral(Q55f,(NL/2-k)*h/NL,(NL/2-k+1)*h/NL,'ArrayValued',true);

A11=A11+a11;
A22=A22+a22;
A12=A12+a12;
A21=A21+a21;
A66=A66+a66;
B11=B11+b11;
B22=B22+b22;
B12=B12+b12;
B21=B21+b21;
B66=B66+b66;
D11=D11+d11;
D22=D22+d22;
D12=D12+d12;
D21=D21+d21;
D66=D66+d66;
E11=E11+e11;
E22=E22+e22;
E12=E12+e12;
E21=E21+e21;
E66=E66+e66;
F11=F11+f11;
F22=F22+f22;
F12=F12+f12;
F21=F21+f21;
F66=F66+f66;
H11=H11+h11;
H22=H22+h22;
H12=H12+h12;
H21=H21+h21;
H66=H66+h66;
A44=A44+a44;
A55=A55+a55;
D44=D44+d44;
D55=D55+d55;
F44=F44+f44;
F55=F55+f55;

NxxT=NxxT+nxt; NyyT=NyyT+nyt; MxxT=MxxT+mxt; MyyT=MyyT+myt; PxxT=PxxT+pxt; PyyT=PyyT+pyt;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
