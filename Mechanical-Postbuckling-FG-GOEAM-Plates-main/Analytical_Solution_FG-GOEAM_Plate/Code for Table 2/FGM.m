%% Material properties
T=300;
T0=300;
DeltaT=T-T0;

% Metal (Al) - Ceramic (Al2O3)
Ec =380E9;
Nuc = 0.3;
Alphac =7.4e-6;

Em =70E9;
Num =0.3;
Alpham =23e-6;

%% FGM functions
E = @(z) Em + (Ec - Em) .* (h.^(-1) .* ((1./2) * h + z)).^p;
Nu = @(z) Num + (Nuc - Num) .* (h.^(-1) .* ((1./2) .* h + z)).^p;
Alpha = @(z) Alpham + (Alphac - Alpham) .* (h.^(-1) .* ((1./2) .* h + z)).^p;

%% Qij matrix
Q11a=@(z) (E(z)./(1-Nu(z).^2)); Q22a=@(z) (E(z)./(1-Nu(z).^2)); Q12a=@(z) (E(z).*Nu(z)./(1-Nu(z).^2));
Q21a=@(z) (E(z).*Nu(z)./(1-Nu(z).^2)); Q66a=@(z) (E(z)./2./(1+Nu(z)));

Q11b=@(z) (E(z)./(1-Nu(z).^2)).*z; Q22b=@(z) (E(z)./(1-Nu(z).^2)).*z; Q12b=@(z) (E(z).*Nu(z)./(1-Nu(z).^2)).*z;
Q21b=@(z) (E(z).*Nu(z)./(1-Nu(z).^2)).*z; Q66b=@(z) (E(z)./2./(1+Nu(z))).*z;

Q11d=@(z) (E(z)./(1-Nu(z).^2)).*(z.^2); Q22d=@(z) (E(z)./(1-Nu(z).^2)).*(z.^2); Q12d=@(z) (E(z).*Nu(z)./(1-Nu(z).^2)).*(z.^2);
Q21d=@(z) (E(z).*Nu(z)./(1-Nu(z).^2)).*(z.^2); Q66d=@(z) (E(z)./2./(1+Nu(z))).*(z.^2);

Q11e=@(z) (E(z)./(1-Nu(z).^2)).*(z.^3); Q22e=@(z) (E(z)./(1-Nu(z).^2)).*(z.^3); Q12e=@(z) (E(z).*Nu(z)./(1-Nu(z).^2)).*(z.^3);
Q21e=@(z) (E(z).*Nu(z)./(1-Nu(z).^2)).*(z.^3); Q66e=@(z) (E(z)./2./(1+Nu(z))).*(z.^3);

Q11f=@(z) (E(z)./(1-Nu(z).^2)).*(z.^4); Q22f=@(z) (E(z)./(1-Nu(z).^2)).*(z.^4); Q12f=@(z) (E(z).*Nu(z)./(1-Nu(z).^2)).*(z.^4);
Q21f=@(z) (E(z).*Nu(z)./(1-Nu(z).^2)).*(z.^4); Q66f=@(z) (E(z)./2./(1+Nu(z))).*(z.^4);

Q11h=@(z) (E(z)./(1-Nu(z).^2)).*(z.^6); Q22h=@(z) (E(z)./(1-Nu(z).^2)).*(z.^6); Q12h=@(z) (E(z).*Nu(z)./(1-Nu(z).^2)).*(z.^6);
Q21h=@(z) (E(z).*Nu(z)./(1-Nu(z).^2)).*(z.^6); Q66h=@(z) (E(z)./2./(1+Nu(z))).*(z.^6);

Q44a=@(z) (E(z)./2./(1+Nu(z))); Q55a=@(z) (E(z)./2./(1+Nu(z)));
Q44d=@(z) (E(z)./2./(1+Nu(z))).*(z.^2); Q55d=@(z) (E(z)./2./(1+Nu(z))).*(z.^2);
Q44f=@(z) (E(z)./2./(1+Nu(z))).*(z.^4); Q55f=@(z) (E(z)./2./(1+Nu(z))).*(z.^4);

nx=@(z) (Q11a(z)+Q12a(z)).*Alpha(z).*DeltaT; ny=@(z) (Q21a(z)+Q22a(z)).*Alpha(z).*DeltaT;
mx=@(z) (Q11a(z)+Q12a(z)).*Alpha(z).*DeltaT.*z; my=@(z) (Q21a(z)+Q22a(z)).*Alpha(z).*DeltaT.*z;
px=@(z) (Q11a(z)+Q12a(z)).*Alpha(z).*DeltaT.*(z.^3); py=@(z) (Q21a(z)+Q22a(z)).*Alpha(z).*DeltaT.*(z.^3);

%%%%%%%%%%%%%%%%%  plate stiffness coefficients %%%%%%%%%%%%%%%%%%%%%%%%%%%
NxxT=integral(nx, -h/2, h/2); NyyT=integral(ny, -h/2, h/2);
MxxT=integral(mx, -h/2, h/2); MyyT=integral(my, -h/2, h/2);
PxxT=integral(px, -h/2, h/2); PyyT=integral(py, -h/2, h/2);

A11=integral(Q11a, -h/2, h/2); A12=integral(Q12a, -h/2, h/2);
A21=integral(Q21a, -h/2, h/2); A22=integral(Q22a, -h/2, h/2);
A66=integral(Q66a, -h/2, h/2);

B11=integral(Q11b, -h/2, h/2); B12=integral(Q12b, -h/2, h/2);
B21=integral(Q21b, -h/2, h/2); B22=integral(Q22b, -h/2, h/2);
B66=integral(Q66b, -h/2, h/2);

D11=integral(Q11d, -h/2, h/2); D12=integral(Q12d, -h/2, h/2);
D21=integral(Q21d, -h/2, h/2); D22=integral(Q22d, -h/2, h/2);
D66=integral(Q66d, -h/2, h/2);

E11=integral(Q11e, -h/2, h/2); E12=integral(Q12e, -h/2, h/2);
E21=integral(Q21e, -h/2, h/2); E22=integral(Q22e, -h/2, h/2);
E66=integral(Q66e, -h/2, h/2);

F11=integral(Q11f, -h/2, h/2); F12=integral(Q12f, -h/2, h/2);
F21=integral(Q21f, -h/2, h/2); F22=integral(Q22f, -h/2, h/2);
F66=integral(Q66f, -h/2, h/2);

H11=integral(Q11h, -h/2, h/2); H12=integral(Q12h, -h/2, h/2);
H21=integral(Q21h, -h/2, h/2); H22=integral(Q22h, -h/2, h/2);
H66=integral(Q66h, -h/2, h/2);

A44=integral(Q44a, -h/2, h/2); A55=integral(Q55a, -h/2, h/2);
D44=integral(Q44d, -h/2, h/2); D55=integral(Q55d, -h/2, h/2);
F44=integral(Q44f, -h/2, h/2); F55=integral(Q55f, -h/2, h/2);
