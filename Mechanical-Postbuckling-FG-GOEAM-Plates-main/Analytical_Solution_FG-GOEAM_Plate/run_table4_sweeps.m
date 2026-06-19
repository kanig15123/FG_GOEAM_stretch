clear all; clc;

% Geometrical parameters
h = 0.005; 
var1 = 10; % a/h
var2 = 1.0; % b/a
a = var1*h;
b = var2*a;
m = 1; n = 1;
c1 = 4/h^2/3; c2 = 3*c1;

% Foundation is zero
Kun = 0; Ksn = 0; Kln = 0;
Ku=0; Kl=0; Ks=0; K1=0; K2=0;

% Loading: Biaxial Compression
Chi1 = -1; Chi2 = -1;

% Table 4 Part 1: HGr = 80%, vary WGr = [0, 0.5, 1.0, 1.5, 2.0, 2.5] wt%
WGr_list = [0, 0.5, 1.0, 1.5, 2.0, 2.5] / 100;
HGr = 80/100;
T = 300;

fprintf('========================================================================\n');
fprintf('TABLE 4 PART 1: HGr = 80%%, VARYING WGr (Biaxial, a/h=10, b/a=1)\n');
fprintf('========================================================================\n');
fprintf('Pattern   WGr(wt%%)   No-Stretch   Stretch       Diff        %% Change\n');
fprintf('------------------------------------------------------------------------\n');

for pattern = 1:3
    if pattern == 1; pat_lbl = 'U-WGr';
    elseif pattern == 2; pat_lbl = 'X-WGr';
    else; pat_lbl = 'O-WGr';
    end
    
    for WGr = WGr_list
        [Ncr_ns, Ncr_s] = solve_plate_Ncr(h, a, b, m, n, c1, c2, T, HGr, WGr, pattern, Chi1, Chi2);
        diff = Ncr_ns - Ncr_s;
        pct = (diff / Ncr_ns) * 100;
        fprintf('%-7s   %4.1f%%      %.6f     %.6f     %.6f    %.4f%%\n', pat_lbl, WGr*100, Ncr_ns, Ncr_s, diff, pct);
    end
end

% Table 4 Part 2: WGr = 2.5 wt%, vary HGr = [0, 20, 40, 60, 80, 100] %
HGr_list = [0, 20, 40, 60, 80, 100] / 100;
WGr = 2.5/100;

fprintf('\n========================================================================\n');
fprintf('TABLE 4 PART 2: WGr = 2.5 wt%%, VARYING HGr (Biaxial, a/h=10, b/a=1)\n');
fprintf('========================================================================\n');
fprintf('Pattern   HGr(%%)     No-Stretch   Stretch       Diff        %% Change\n');
fprintf('------------------------------------------------------------------------\n');

for pattern = 1:3
    if pattern == 1; pat_lbl = 'U-WGr';
    elseif pattern == 2; pat_lbl = 'X-WGr';
    else; pat_lbl = 'O-WGr';
    end
    
    for HGr = HGr_list
        [Ncr_ns, Ncr_s] = solve_plate_Ncr(h, a, b, m, n, c1, c2, T, HGr, WGr, pattern, Chi1, Chi2);
        diff = Ncr_ns - Ncr_s;
        pct = (diff / Ncr_ns) * 100;
        fprintf('%-7s   %4.0f%%      %.6f     %.6f     %.6f    %.4f%%\n', pat_lbl, HGr*100, Ncr_ns, Ncr_s, diff, pct);
    end
end


function [Ncr_ns, Ncr_s] = solve_plate_Ncr(h, a, b, m, n, c1, c2, T, HGr, WGr, pattern, Chi1, Chi2)
    % Material Properties
    T0 = 300; DeltaT = T - T0; 
    EGr = 929.57e9; NuGr = 0.22; AlphaGr=-3.98e-6; RoGr=1800;
    ECu = 65.79e9; NuCu = 0.387; AlphaCu=16.51e-6; RoCu=8800;
    Dm = ECu * h^3 / 12 / (1 - NuCu^2);
    NL = 10;
    VGrs = WGr / (WGr + (RoGr/RoCu)*(1 - WGr));
    
    %% --- NO STRETCH ---
    A11=0; A22=0; A12=0; A21=0; A66=0;
    B11=0; B22=0; B12=0; B21=0; B66=0;
    D11=0; D22=0; D12=0; D21=0; D66=0;
    E11=0; E22=0; E12=0; E21=0; E66=0;
    F11=0; F22=0; F12=0; F21=0; F66=0;
    H11=0; H22=0; H12=0; H21=0; H66=0;
    A44=0; A55=0; D44=0; D55=0; F44=0; F55=0;
    
    for k = 1:NL
        z_bot = (NL/2-k)*h/NL; z_top = (NL/2-k+1)*h/NL;
        z_mid = (z_bot + z_top)/2;
        
        if pattern == 1 % U-WGr
            VGr = VGrs;
        elseif pattern == 2 % X-WGr
            VGr = 2 * VGrs * (2 * abs(z_mid) / h);
        else % O-WGr
            VGr = 2 * VGrs * (1 - 2 * abs(z_mid) / h);
        end
        
        TR = T / T0;
        fE = 1.11 - 1.22*VGr - 0.134*TR + 0.559*VGr*TR - 5.5*HGr*VGr + 38*HGr*VGr^2 - 20.6*HGr^2*VGr^2;
        fv = 1.01 - 1.43*VGr + 0.165*TR - 16.8*HGr*VGr - 1.1*HGr*VGr*TR + 16*HGr^2*VGr^2;
        fAlpha = 0.794 - 16.8*VGr^2 - 0.0279*TR^2 + 0.182*TR*(1+VGr);
        
        lGr=83.76e-10; tGr=3.4e-10; Xi=2*lGr/tGr; Eta=((EGr/ECu)-1)/((EGr/ECu)+Xi); VCu=1-VGr;
        E_l = (1+Xi*Eta*VGr)*ECu*fE/(1-Eta*VGr); Nu_l = (NuGr*VGr+NuCu*VCu)*fv; Alpha_l = (AlphaGr*VGr+AlphaCu*VCu)*fAlpha;
        
        constQ = 1 - Nu_l^2;
        Q11_val = E_l / constQ; Q12_val = E_l * Nu_l / constQ; Q66_val = E_l / (2*(1+Nu_l)); Q44_val = Q66_val;
        
        % Integrals
        a11 = Q11_val * (z_top - z_bot);
        a12 = Q12_val * (z_top - z_bot);
        a66 = Q66_val * (z_top - z_bot);
        
        b11 = Q11_val * (z_top^2 - z_bot^2)/2;
        b12 = Q12_val * (z_top^2 - z_bot^2)/2;
        b66 = Q66_val * (z_top^2 - z_bot^2)/2;
        
        d11 = Q11_val * (z_top^3 - z_bot^3)/3;
        d12 = Q12_val * (z_top^3 - z_bot^3)/3;
        d66 = Q66_val * (z_top^3 - z_bot^3)/3;
        
        e11 = Q11_val * (z_top^4 - z_bot^4)/4;
        e12 = Q12_val * (z_top^4 - z_bot^4)/4;
        e66 = Q66_val * (z_top^4 - z_bot^4)/4;
        
        f11 = Q11_val * (z_top^5 - z_bot^5)/5;
        f12 = Q12_val * (z_top^5 - z_bot^5)/5;
        f66 = Q66_val * (z_top^5 - z_bot^5)/5;
        
        h11 = Q11_val * (z_top^7 - z_bot^7)/7;
        h12 = Q12_val * (z_top^7 - z_bot^7)/7;
        h66 = Q66_val * (z_top^7 - z_bot^7)/7;
        
        a44 = Q44_val * (z_top - z_bot);
        d44 = Q44_val * (z_top^3 - z_bot^3)/3;
        f44 = Q44_val * (z_top^5 - z_bot^5)/5;
        
        A11=A11+a11; A22=A22+a11; A12=A12+a12; A21=A21+a12; A66=A66+a66;
        B11=B11+b11; B22=B22+b11; B12=B12+b12; B21=B21+b12; B66=B66+b66;
        D11=D11+d11; D22=D22+d11; D12=D12+d12; D21=D21+d12; D66=D66+d66;
        E11=E11+e11; E22=E22+e11; E12=E12+e12; E21=E21+e12; E66=E66+e66;
        F11=F11+f11; F22=F22+f11; F12=F12+f12; F21=F21+f12; F66=F66+f66;
        H11=H11+h11; H22=H22+h11; H12=H12+h12; H21=H21+h12; H66=H66+h66;
        A44=A44+a44; A55=A55+a44; D44=D44+d44; D55=D55+d44; F44=F44+f44; F55=F55+f44;
    end
    
    S2 = 0; % No foundation
    
    % Solver for no-stretch
    lamda0 = solve_lamda0(A11, A22, A12, A21, A66, B11, B22, B12, B21, B66, ...
                         D11, D22, D12, D21, D66, E11, E22, E12, E21, E66, ...
                         F11, F22, F12, F21, F66, H11, H22, H12, H21, H66, ...
                         A44, A55, D44, D55, F44, F55, S2, a, b, m, n, c1, c2, Chi1, Chi2);
                     
    Ncr_ns = lamda0 * a^2 / (ECu * h^3);
    
    %% --- WITH STRETCH (STATIC CONDENSATION) ---
    A33=0; A13=0; A23=0;
    
    for k = 1:NL
        z_bot = (NL/2-k)*h/NL; z_top = (NL/2-k+1)*h/NL;
        z_mid = (z_bot + z_top)/2;
        
        if pattern == 1 % U-WGr
            VGr = VGrs;
        elseif pattern == 2 % X-WGr
            VGr = 2 * VGrs * (2 * abs(z_mid) / h);
        else % O-WGr
            VGr = 2 * VGrs * (1 - 2 * abs(z_mid) / h);
        end
        
        TR = T / T0;
        fE = 1.11 - 1.22*VGr - 0.134*TR + 0.559*VGr*TR - 5.5*HGr*VGr + 38*HGr*VGr^2 - 20.6*HGr^2*VGr^2;
        fv = 1.01 - 1.43*VGr + 0.165*TR - 16.8*HGr*VGr - 1.1*HGr*VGr*TR + 16*HGr^2*VGr^2;
        
        lGr=83.76e-10; tGr=3.4e-10; Xi=2*lGr/tGr; Eta=((EGr/ECu)-1)/((EGr/ECu)+Xi); VCu=1-VGr;
        E_l = (1+Xi*Eta*VGr)*ECu*fE/(1-Eta*VGr); Nu_l = (NuGr*VGr+NuCu*VCu)*fv;
        
        constC = (1 + Nu_l) * (1 - 2*Nu_l);
        C12_3D_val = E_l * Nu_l / constC; C33_3D_val = E_l * (1 - Nu_l) / constC;
        
        a33 = C33_3D_val * (z_top - z_bot);
        a13 = C12_3D_val * (z_top - z_bot);
        
        A33 = A33 + a33;
        A13 = A13 + a13;
    end
    A23 = A13;
    
    % Static condensation on plane stiffnesses
    A11_cond = A11 - (A13^2)/A33;
    A22_cond = A22 - (A23^2)/A33;
    A12_cond = A12 - (A13*A23)/A33;
    A21_cond = A12_cond;
    
    lamda0_s = solve_lamda0(A11_cond, A22_cond, A12_cond, A21_cond, A66, B11, B22, B12, B21, B66, ...
                           D11, D22, D12, D21, D66, E11, E22, E12, E21, E66, ...
                           F11, F22, F12, F21, F66, H11, H22, H12, H21, H66, ...
                           A44, A55, D44, D55, F44, F55, S2, a, b, m, n, c1, c2, Chi1, Chi2);
                       
    Ncr_s = lamda0_s * a^2 / (ECu * h^3);
end

function lamda0 = solve_lamda0(A11, A22, A12, A21, A66, B11, B22, B12, B21, B66, ...
                              D11, D22, D12, D21, D66, E11, E22, E12, E21, E66, ...
                              F11, F22, F12, F21, F66, H11, H22, H12, H21, H66, ...
                              A44, A55, D44, D55, F44, F55, S2, a, b, m, n, c1, c2, Chi1, Chi2)
    Denom = (A12.*A21 - A11.*A22);
    C11 = A12 ./ Denom;
    C12 = -(A22 ./ Denom);
    C13 = -((-A22.*B11 + A12.*B21) ./ Denom);
    C14 = -((-A22.*B12 + A12.*B22) ./ Denom);
    C15 = -((-A22.*E11 + A12.*E21) ./ Denom);
    C16 = -((-A22.*E12 + A12.*E22) ./ Denom);

    C21 = A11 ./ (-Denom);
    C22 = -(A21 ./ (-Denom));
    C23 = -((-A21.*B11 + A11.*B21) ./ (-Denom));
    C24 = -((-A21.*B12 + A11.*B22) ./ (-Denom));
    C25 = -((-A21.*E11 + A11.*E21) ./ (-Denom));
    C26 = -((-A21.*E12 + A11.*E22) ./ (-Denom));

    C31 = -(1./A66);
    C32 = -(B66./A66);
    C33 = -(E66./A66);

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

    Denominator_R = ((b^4 * J11 * m^4 + a^2 * b^2 * J12 * m^2 * n^2 + a^4 * J13 * n^4) * pi);

    R11 = (a^2 * n^2) / (32 * b^2 * J11 * m^2);
    R21 = (b^2 * m^2) / (32 * a^2 * J13 * n^2);
    R31 = (-a * b^4 * J21 * m^3 - a^3 * b^2 * J22 * m * n^2) / Denominator_R;
    R32 = (-a^2 * b^3 * J23 * m^2 * n - a^4 * b * J24 * n^3) / Denominator_R;
    R33 = (-b^4 * J31 * m^4 * pi - a^2 * b^2 * J32 * m^2 * n^2 * pi - a^4 * J33 * n^4 * pi) / Denominator_R;

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

    Denom_P = (V24 * V33 - V23 * V34);
    P11 = -((V24 * V31 - V21 * V34) / Denom_P);
    P21 = -((-V23 * V31 + V21 * V33) / Denom_P);

    Term_Common = (Chi1.*V18 + Chi2.*V19).^(-1);
    lamda0 = (-1) .* (S2 + V11 + P11.*V14 + P21.*V15) .* Term_Common;
end
