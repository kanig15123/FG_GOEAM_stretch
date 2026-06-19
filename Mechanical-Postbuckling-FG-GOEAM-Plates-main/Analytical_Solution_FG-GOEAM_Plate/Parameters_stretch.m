%% ========================================================================
%  File: Parameters_Analytical_Solution.m (QUASI-3D CORRECTED)
%  Description: Analytical parameters with Thickness Stretching Correction
% =========================================================================

%% ---------------- QUASI-3D (THICKNESS STRETCHING) CORRECTION ----------------
% Bu blok, FG_GOEAM dosyasindan gelen A33 ve A13 terimlerini kullanarak
% "Plane Stress" matrislerini "Quasi-3D" efektif matrislerine donusturur.
% Mantik: Ez (epsilon_z) etkisi, duzlem ici rijitligi azaltir (relaxation).

if exist('A11_3D', 'var') && any(A11_3D(:) ~= 0)
    % 1. Mekanik Rijitlik Indirgemesi (Static Condensation)
    detK = A11_3D .* D11_3D - B11_3D.^2;
    I11 = D11_3D ./ detK;
    I12 = -B11_3D ./ detK;
    I21 = I12;
    I22 = A11_3D ./ detK;

    % A condensed:
    A11 = A11_3D - A12_3D.*(I11.*A12_3D + I12.*B12_3D) - B12_3D.*(I21.*A12_3D + I22.*B12_3D);
    A22 = A11;
    A12 = A12_3D - A12_3D.*(I11.*A12_3D + I12.*B12_3D) - B12_3D.*(I21.*A12_3D + I22.*B12_3D);
    A21 = A12;

    % B condensed:
    B11 = B11_3D - A12_3D.*(I11.*B12_3D + I12.*D12_3D) - B12_3D.*(I21.*B12_3D + I22.*D12_3D);
    B22 = B11;
    B12 = B12_3D - A12_3D.*(I11.*B12_3D + I12.*D12_3D) - B12_3D.*(I21.*B12_3D + I22.*D12_3D);
    B21 = B12;

    % D condensed:
    D11 = D11_3D - B12_3D.*(I11.*B12_3D + I12.*D12_3D) - D12_3D.*(I21.*B12_3D + I22.*D12_3D);
    D22 = D11;
    D12 = D12_3D - B12_3D.*(I11.*B12_3D + I12.*D12_3D) - D12_3D.*(I21.*B12_3D + I22.*D12_3D);
    D21 = D12;

    % E condensed:
    E11 = E11_3D - A12_3D.*(I11.*E12_3D + I12.*F12_3D) - B12_3D.*(I21.*E12_3D + I22.*F12_3D);
    E22 = E11;
    E12 = E12_3D - A12_3D.*(I11.*E12_3D + I12.*F12_3D) - B12_3D.*(I21.*E12_3D + I22.*F12_3D);
    E21 = E12;

    % F condensed:
    F11 = F11_3D - B12_3D.*(I11.*E12_3D + I12.*F12_3D) - D12_3D.*(I21.*E12_3D + I22.*F12_3D);
    F22 = F11;
    F12 = F12_3D - B12_3D.*(I11.*E12_3D + I12.*F12_3D) - D12_3D.*(I21.*E12_3D + I22.*F12_3D);
    F21 = F12;

    % H condensed:
    H11 = H11_3D - E12_3D.*(I11.*E12_3D + I12.*F12_3D) - F12_3D.*(I21.*E12_3D + I22.*F12_3D);
    H22 = H11;
    H12 = H12_3D - E12_3D.*(I11.*E12_3D + I12.*F12_3D) - F12_3D.*(I21.*E12_3D + I22.*F12_3D);
    H21 = H12;

    % 2. Termal Kuvvet Duzeltmesi
    if exist('NzzT', 'var')
        NxxT = NxxT - A12_3D.*(I11.*NzzT) - B12_3D.*(I22.*NzzT);
        NyyT = NyyT - A12_3D.*(I11.*NzzT) - B12_3D.*(I22.*NzzT);
    end
else
    % A33 yoksa standart HSDT (Plane Stress) devam eder.
    % disp('Standard HSDT (No Thickness Stretching) used.');
end

%% Eq. (19) - Inverted Stiffness Coefficients
% Yukarida A matrisleri guncellendigi icin bu denklemler artik
% otomatik olarak Quasi-3D etkisini icerir. Formulleri degistirmeye gerek yoktur.

Denom = (A12.*A21 - A11.*A22); % Payda (Determinant benzeri)

C11 = A12 ./ Denom;
C12 = -(A22 ./ Denom);
C13 = -((-A22.*B11 + A12.*B21) ./ Denom);
C14 = -((-A22.*B12 + A12.*B22) ./ Denom);
C15 = -((-A22.*E11 + A12.*E21) ./ Denom);
C16 = -((-A22.*E12 + A12.*E22) ./ Denom);

C21 = A11 ./ (-Denom); % Dikkat: Isaret farki (A11*A22 - A12*A21)
C22 = -(A21 ./ (-Denom));
C23 = -((-A21.*B11 + A11.*B21) ./ (-Denom));
C24 = -((-A21.*B12 + A11.*B22) ./ (-Denom));
C25 = -((-A21.*E11 + A11.*E21) ./ (-Denom));
C26 = -((-A21.*E12 + A11.*E22) ./ (-Denom));

% Shear & Higher Order Terms
C31 = -(1./A66);
C32 = -(B66./A66);
C33 = -(E66./A66);

%% Eq. (21) - Geometric Coefficients
J11 = C21;
J12 = C11 + C22 - C31;
J13 = C12;
J21 = C23 - c1 * C25;
J22 = C13 - c1 * C15 - C32 + c1 * C33;
J23 = C24 - c1 * C26 - C32 + c1 * C33;
J24 = C14 - c1 * C16;
J31 = -c1 * C25;
J32 = -c1 * C15 - c1 * C26 + 2 * c1 * C33;
J33 = -c1 * C16;

%% Appendix A - Governing Equation Coefficients
S11 = -c1^2*C15*E11 - c1^2*C25*E12 - c1^2*H11;
S12 = -c1^2*C16*E11 - c1^2*C26*E12 - c1^2*C15*E21 - c1^2*C25*E22 - 4*c1^2*C33*E66 - c1^2*H12 - c1^2*H21 - 4*c1^2*H66;
S13 = -c1^2*C16*E21 - c1^2*C26*E22 - c1^2*H22;
S14 = A55 - 2*c2*D55 + c2^2*F55;
S15 = A44 - 2*c2*D44 + c2^2*F44;
S16 = A55 - 2*c2*D55 + c2^2*F55;
S17 = c1*C13*E11 - c1^2*C15*E11 + c1*C23*E12 - c1^2*C25*E12 + c1*F11 - c1^2*H11;
S18 = c1*C13*E21 - c1^2*C15*E21 + c1*C23*E22 - c1^2*C25*E22 + 2*c1*C32*E66 - 2*c1^2*C33*E66 + c1*F21 + 2*c1*F66 - c1^2*H21 - 2*c1^2*H66;
S19 = A44 - 2*c2*D44 + c2^2*F44;
S110 = c1*C14*E21 - c1^2*C16*E21 + c1*C24*E22 - c1^2*C26*E22 + c1*F22 - c1^2*H22;
S111 = c1*C14*E11 - c1^2*C16*E11 + c1*C24*E12 - c1^2*C26*E12 + 2*c1*C32*E66 - 2*c1^2*C33*E66 + c1*F12 + 2*c1*F66 - c1^2*H12 - 2*c1^2*H66;
S112 = c1*C11*E11 + c1*C21*E12;
S113 = c1*C12*E11 + c1*C22*E12 + c1*C11*E21 + c1*C21*E22 + 2*c1*C31*E66;
S114 = c1*C12*E21 + c1*C22*E22;

S21 = -B11 * c1 * C16 - B12 * c1 * C26 - 2 * B66 * c1 * C33 + c1^2 * C16 * E11 + ...
      c1^2 * C26 * E12 + 2 * c1^2 * C33 * E66 - c1 * F12 - 2 * c1 * F66 + c1^2 * H12 + ...
      2 * c1^2 * H66;
S22 = -B11 * c1 * C15 - B12 * c1 * C25 + c1^2 * C15 * E11 + c1^2 * C25 * E12 - ...
      c1 * F11 + c1^2 * H11;
S23 = -A55 + 2 * c2 * D55 - c2^2 * F55;
S24 = B11 * C13 - B11 * c1 * C15 + B12 * C23 - B12 * c1 * C25 + D11 - c1 * C13 * E11 + ...
      c1^2 * C15 * E11 - c1 * C23 * E12 + c1^2 * C25 * E12 - 2 * c1 * F11 + c1^2 * H11;
S25 = B66 * C32 - B66 * c1 * C33 + D66 - c1 * C32 * E66 + c1^2 * C33 * E66 - ...
      2 * c1 * F66 + c1^2 * H66;
S26 = -A55 + 2 * c2 * D55 - c2^2 * F55;
S27 = B11 * C14 - B11 * c1 * C16 + B12 * C24 - B12 * c1 * C26 + B66 * C32 - ...
      B66 * c1 * C33 + D12 + D66 - c1 * C14 * E11 + c1^2 * C16 * E11 - c1 * C24 * E12 + ...
      c1^2 * C26 * E12 - c1 * C32 * E66 + c1^2 * C33 * E66 - 2 * c1 * F12 - 2 * c1 * F66 + ...
      c1^2 * H12 + c1^2 * H66;
S28 = B11 * C11 + B12 * C21 - c1 * C11 * E11 - c1 * C21 * E12;
S29 = B11 * C12 + B12 * C22 + B66 * C31 - c1 * C12 * E11 - c1 * C22 * E12 - ...
      c1 * C31 * E66;
      
S31 = -B21 * c1 * C15 - B22 * c1 * C25 - 2 * B66 * c1 * C33 + c1^2 * C15 * E21 + ...
      c1^2 * C25 * E22 + 2 * c1^2 * C33 * E66 - c1 * F21 - 2 * c1 * F66 + c1^2 * H21 + ...
      2 * c1^2 * H66;
S32 = -B21 * c1 * C16 - B22 * c1 * C26 + c1^2 * C16 * E21 + c1^2 * C26 * E22 - ...
      c1 * F22 + c1^2 * H22;
S33 = -A44 + 2 * c2 * D44 - c2^2 * F44;
S34 = B21 * C13 - B21 * c1 * C15 + B22 * C23 - B22 * c1 * C25 + B66 * C32 - ...
      B66 * c1 * C33 + D21 + D66 - c1 * C13 * E21 + c1^2 * C15 * E21 - c1 * C23 * E22 + ...
      c1^2 * C25 * E22 - c1 * C32 * E66 + c1^2 * C33 * E66 - 2 * c1 * F21 - 2 * c1 * F66 + ...
      c1^2 * H21 + c1^2 * H66;
S35 = B66 * C32 - B66 * c1 * C33 + D66 - c1 * C32 * E66 + c1^2 * C33 * E66 - ...
      2 * c1 * F66 + c1^2 * H66;
S36 = B21 * C14 - B21 * c1 * C16 + B22 * C24 - B22 * c1 * C26 + D22 - c1 * C14 * E21 + ...
      c1^2 * C16 * E21 - c1 * C24 * E22 + c1^2 * C26 * E22 - 2 * c1 * F22 + c1^2 * H22;
S37 = -A44 + 2 * c2 * D44 - c2^2 * F44;
S38 = B21 * C12 + B22 * C22 - c1 * C12 * E21 - c1 * C22 * E22;
S39 = B21 * C11 + B22 * C21 + B66 * C31 - c1 * C11 * E21 - c1 * C21 * E22 - ...
      c1 * C31 * E66;

%% Eq. (28) - R Coefficients (Stress Function)
Denominator_R = ((b^4 * J11 * m^4 + a^2 * b^2 * J12 * m^2 * n^2 + a^4 * J13 * n^4) * pi);

R11 = (a^2 * n^2) / (32 * b^2 * J11 * m^2);
R21 = (b^2 * m^2) / (32 * a^2 * J13 * n^2);
R31 = (-a * b^4 * J21 * m^3 - a^3 * b^2 * J22 * m * n^2) / Denominator_R;
R32 = (-a^2 * b^3 * J23 * m^2 * n - a^4 * b * J24 * n^3) / Denominator_R;
R33 = (-b^4 * J31 * m^4 * pi - a^2 * b^2 * J32 * m^2 * n^2 * pi - a^4 * J33 * n^4 * pi) / Denominator_R;

%% Eq. (32) - Foundation (S2)
S2=(1/8).*a.^(-1).*b.^(-1).*m.^(-1).*pi.^(-1).*(b.^2.*K2.*m.^2.*pi.^2+ ...
  a.^2.*(b.^2.*K1+K2.*n.^2.*pi.^2)).*((-2).*m.*pi+sin(2.*m.*pi));

%% Eq. (31) - Galerkin Integrals (V Coefficients)
% Bu integraller geometri ve mod sekillerine bagli oldugu icin degismez.
V11 = 1 / (8 * a^3 * b^3 * m) * pi * (a^4 * n^4 * pi^2 * (R33 * S114 + S13) + ...
    b^4 * (m^4 * pi^2 * (S11 + R33 * S112) - a^2 * m^2 * S14) + ...
    a^2 * b^2 * n^2 * (m^2 * pi^2 * (R33 * S113 + S12) - a^2 * S15)) * (2 * m * pi - sin(2 * m * pi));
V12 = -(1 / (3 * a^3 * b^3 * m * n)) * 4 * pi^2 * (-((1 + 2 * (-1)^n) * a^2 * n^2 * (-b^2 * m^2 * R33 + ...
    16 * a^2 * n^2 * R21 * S114)) + b^2 * m^2 * (a^2 * n^2 * R33 - 16 * b^2 * m^2 * R11 * S112) * ...
    (2 * cos(m * pi) + cos(2 * m * pi))) * sin((m * pi) / 2)^2 * sin((n * pi) / 2)^2;
V13 = -1 / (8 * a * b) * m * n^2 * pi^3 * (4 * m * pi * (R11 + R21) - ...
    2 * (2 * R11 + R21) * sin(2 * m * pi) + R11 * sin(4 * m * pi));
V14 = -1 / (8 * a^3 * b^3 * m) * (a^4 * n^4 * pi^3 * R31 * S114 + ...
    b^4 * (m^4 * pi^3 * R31 * S112 - a^3 * m * S16 + a * m^3 * pi^2 * S17) + ...
    a^2 * b^2 * m * n^2 * pi^2 * (m * pi * R31 * S113 + a * S18)) * (-2 * m * pi + sin(2 * m * pi));
V15 = -1 / (8 * a^3 * b^3 * m) * (pi^2 * (b^4 * m^4 * pi * R32 * S112 + ...
    a^2 * b^2 * m^2 * n * (b * S111 + n * pi * R32 * S113) + ...
    a^4 * n^3 * (b * S110 + n * pi * R32 * S114)) - a^4 * b^3 * n * S19) * (-2 * m * pi + sin(2 * m * pi));
V16 = -1 / (3 * a * b) * 4 * m * n * pi^2 * R31 * (1 + 2 * (-1)^n + 2 * cos(m * pi) + cos(2 * m * pi)) * ...
    sin((m * pi) / 2)^2 * sin((n * pi) / 2)^2;
V17 = -1 / (3 * a * b) * 4 * m * n * pi^2 * R32 * (1 + 2 * (-1)^n + 2 * cos(m * pi) + cos(2 * m * pi)) * ...
    sin((m * pi) / 2)^2 * sin((n * pi) / 2)^2;
V18 = -(b * m * pi * (2 * m * pi - sin(2 * m * pi))) / (8 * a);
V19 = -(a * n^2 * pi * (2 * m * pi - sin(2 * m * pi))) / (8 * b * m);
V21 = -1 / (8 * a^2 * b) * (b^2 * m^2 * pi^2 * (S22 + R33 * S28) + ...
    a^2 * (-b^2 * S23 + n^2 * pi^2 * (S21 + R33 * S29))) * (2 * m * pi + sin(2 * m * pi));
V22 = (16 * (-1 + (-1)^n) * b * m^2 * pi * R11 * S28 * (-1 + cos(m * pi)^3)) / (3 * a^2 * n);
V23 = -1 / (8 * a^2 * b * m * pi) * (a * b^2 * m^2 * pi^2 * S24 + ...
    a^3 * (n^2 * pi^2 * S25 - b^2 * S26) + b^2 * m^3 * pi^3 * R31 * S28 + ...
    a^2 * m * n^2 * pi^3 * R31 * S29) * (2 * m * pi + sin(2 * m * pi));
V24 = -1 / (8 * a^2 * b) * pi * (b^2 * m^2 * pi * R32 * S28 + ...
    a^2 * n * (b * S27 + n * pi * R32 * S29)) * (2 * m * pi + sin(2 * m * pi));
V31 = 1 / (8 * a * b^2 * m) * n * (a^2 * n^2 * pi^2 * (S32 + R33 * S38) + ...
    b^2 * (-a^2 * S33 + m^2 * pi^2 * (S31 + R33 * S39))) * (-2 * m * pi + sin(2 * m * pi));
V32 = (16 * (-1 + (-1)^n) * a * n^2 * pi * R21 * S38 * (-1 + cos(m * pi))) / (3 * b^2 * m);
V33 = -1 / (8 * a * b^2 * m) * n * pi * (a * b^2 * m * S34 + ...
    a^2 * n^2 * pi * R31 * S38 + b^2 * m^2 * pi * R31 * S39) * (2 * m * pi - sin(2 * m * pi));
V34 = 1 / (8 * a * b^2 * m * pi) * (a^2 * b * n^2 * pi^2 * S36 + ...
    b^3 * (m^2 * pi^2 * S35 - a^2 * S37) + a^2 * n^3 * pi^3 * R32 * S38 + ...
    b^2 * m^2 * n * pi^3 * R32 * S39) * (-2 * m * pi + sin(2 * m * pi));

%% Eq. (34) - Postbuckling Coefficients
Denom_P = (V24 * V33 - V23 * V34);
P11 = -((V24 * V31 - V21 * V34) / Denom_P);
P12 = -((V24 * V32 - V22 * V34) / Denom_P);
P21 = -((-V23 * V31 + V21 * V33) / Denom_P);
P22 = -((-V23 * V32 + V22 * V33) / Denom_P);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% FINAL STABILITY EQUATIONS (Eq. 38-39)
% Lamda0: Kritik Burkulma Yuku (Non-dimensional)
% Lamda2: Postbuckling Egimi

Term_Common = (Chi1.*V18 + Chi2.*V19).^(-1);

lamda0 = (-1) .* (S2 + V11 + P11.*V14 + P21.*V15) .* Term_Common;
lamda1 = (-1) .* h .* (V12 + P12.*V14 + P22.*V15 + P11.*V16 + P21.*V17) .* Term_Common;
lamda2 = (-1) .* h.^2 .* (V13 + P12.*V16 + P22.*V17) .* Term_Common;

% fprintf('Buckling Load (Lambda0): %.4f\n', lamda0);