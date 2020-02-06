% Generates TF pair potential energy tables from 'Startpoint' up to a given
% length input 'Endpoint'. Plotswitch is a logical true or false to
% determine whether to plot the resulting potentials
% Recommended spacing is 0.005 angstroms or 0.0005 nm
% INPUT parameters units are: Energy = kJ/mol, length = nm
% OUTPUT units are: Energy = kJ/mol, length = nm
% These units are used in GROMACS by default.
% CRDamping is a boolean which adds a close range damping to the
% attractive R^-6 term

% C6Damping: Places medium range damping on C6/r^6 terms
% 0 = no (default) damping
% 1 = BJ/rational damping (same as in D3(BJ), damps to a constant. Fairly
% weak damping)
% 2 = Tang Damping (Mid strength damping, damps to zero at r=0)
% 3 = MMDRE Damping function (very weak damping)
% 4 = PAMoC Damping function (weak damping)
% 5 = EHFSK Damping function (strong damping)
% 6 = WY damping function (strongest damping)

function [U_PM_out, U_PP_out, U_MM_out] = TF_Potential_Generator(Startpoint,...
    Endpoint,Spacing,Salt,Parameters,plotswitch,vdw_modifier,RVDW_Cutoff,CRDamping,C6Damping)

%% Conversion factors and fundamental constants
nm_per_m = 1e+9; % nm per m
NA = 6.0221409e23; % Molecules per mole
e_c = 1.60217662e-19; % Elementary charge in Coulombs
epsilon_0 = (8.854187817620e-12)*1000/(nm_per_m*NA); % Vacuum Permittivity C^2 mol kJ^-1 nm^-1
k_0 = 1/(4*pi*epsilon_0); % Coulomb constant in kJ nm C^-2 mol^-1

%% Close-Range Sigmoid Damping Parameters
r_d = 0.10; % If close-range dispersion damping is on, this is the value of the sigmoid's midpoint
b = 150; % sigmoid "steepness" for damping

%% TF Parameter: q (charge)
q.Li =  1; % atomic
q.Na =  1; % atomic
q.K  =  1; % atomic
q.Rb =  1; % atomic
q.Cs =  1; % atomic

q.F  = -1; % atomic
q.Cl = -1; % atomic
q.Br = -1; % atomic
q.I  = -1; % atomic

%% Generate range (r) in nm
r = Startpoint:Spacing:Endpoint;

%% Split Salt Into Component Metal and Halide
[Metal,Halide] = Separate_Metal_Halide(Salt);

%% Parameters from input
alpha.MM = Parameters(1,1);
alpha.XX = Parameters(1,2);
alpha.MX = Parameters(1,3);

B.MM = Parameters(2,1);
B.XX = Parameters(2,2);
B.MX = Parameters(2,3);

C.MM = Parameters(3,1);
C.XX = Parameters(3,2);
C.MX = Parameters(3,3);

D.MM = Parameters(4,1);
D.XX = Parameters(4,2);
D.MX = Parameters(4,3);

%% If Damping at close range
if CRDamping
    f_r = 1./(1 + exp(-b.*(r - r_d))); % sigmoid damping function
    df_r = (b.*exp(-b.*(r - r_d)))./((1 + exp(-b.*(r - r_d))).^2); % sigmoid damping function derivative
    f_cutoff = 1/(1 + exp(-b*(RVDW_Cutoff - r_d))); % damping function value at vdw cutoff
else
    f_r = 1; % No damping
    df_r = 0; % Damping derivative is zero
    f_cutoff = 1; % damping function value at vdw cutoff
end

%% Medium Range Dispersion Damping: No Damping (default)
if C6Damping == 0
    % No-damping function
    f6_MX = f_r;
    f6_MM = f_r;
    f6_XX = f_r;
    
    f8_MX = f_r;
    f8_MM = f_r;
    f8_XX = f_r;

    % Derivative parts of no-damping function   
    df6_MX = df_r;
    df6_MM = df_r;
    df6_XX = df_r;
    
    df8_MX = df_r;
    df8_MM = df_r;
    df8_XX = df_r;

    % Values at vdw cutoff
    f6_cutoff_MX = f_cutoff;
    f6_cutoff_MM = f_cutoff;
    f6_cutoff_XX = f_cutoff;
    f8_cutoff_MX = f_cutoff;
    f8_cutoff_MM = f_cutoff;
    f8_cutoff_XX = f_cutoff;
    
%% BJ-Type Rational Damping
elseif C6Damping == 1

    % Damping distances
    R0_MX = sqrt(D.MX/C.MX); % in nm
    R0_MM = sqrt(D.MM/C.MM); % in nm
    R0_XX = sqrt(D.XX/C.XX); % in nm
    
    % Damping functions (unitless)
    f6_MX = f_r./( 1 + ( R0_MX ./ r ).^6 );
    f6_MM = f_r./( 1 + ( R0_MM ./ r ).^6 );
    f6_XX = f_r./( 1 + ( R0_XX ./ r ).^6 );
    
    f8_MX = f_r./( 1 + ( R0_MX ./ r ).^8 );
    f8_MM = f_r./( 1 + ( R0_MM ./ r ).^8 );
    f8_XX = f_r./( 1 + ( R0_XX ./ r ).^8 );
    
    % Values of damping function at vdw cutoff
    f6_cutoff_MX = f_cutoff/( 1 + ( R0_MX / RVDW_Cutoff )^6 );
    f6_cutoff_MM = f_cutoff/( 1 + ( R0_MM / RVDW_Cutoff )^6 );
    f6_cutoff_XX = f_cutoff/( 1 + ( R0_XX / RVDW_Cutoff )^6 );
    f8_cutoff_MX = f_cutoff/( 1 + ( R0_MX / RVDW_Cutoff )^8 );
    f8_cutoff_MM = f_cutoff/( 1 + ( R0_MM / RVDW_Cutoff )^8 );
    f8_cutoff_XX = f_cutoff/( 1 + ( R0_XX / RVDW_Cutoff )^8 );
    
    % Derivative of damping functions
    df6_MX = f_r.*6.*(R0_MX.^6).*(r.^5)./(((r.^6) + (R0_MX.^6)).^2) + df_r./( 1 + ( R0_MX ./ r ).^6 );
    df6_MM = f_r.*6.*(R0_MM.^6).*(r.^5)./(((r.^6) + (R0_MM.^6)).^2) + df_r./( 1 + ( R0_MM ./ r ).^6 );
    df6_XX = f_r.*6.*(R0_XX.^6).*(r.^5)./(((r.^6) + (R0_XX.^6)).^2) + df_r./( 1 + ( R0_XX ./ r ).^6 );
    
    df8_MX = f_r.*8.*(R0_MX.^8).*(r.^7)./(((r.^8) + (R0_MX.^8)).^2) + df_r./( 1 + ( R0_MX ./ r ).^8 );
    df8_MM = f_r.*8.*(R0_MM.^8).*(r.^7)./(((r.^8) + (R0_MM.^8)).^2) + df_r./( 1 + ( R0_MM ./ r ).^8 );
    df8_XX = f_r.*8.*(R0_XX.^8).*(r.^7)./(((r.^8) + (R0_XX.^8)).^2) + df_r./( 1 + ( R0_XX ./ r ).^8 );
    
%% Tang and Toennies Damping function. Cite:
% "An improved simple model for the van der Waals potential based on universal damping functions for the dispersion coefficients."
% K. T. Tang, J. P. Toennies
% J. Chem. Phys. 1984, 80, 3726-3741.
elseif C6Damping == 2

    % C6 damping functions
    f6sum_MX = 0;
    f6sum_MM = 0;
    f6sum_XX = 0;
    for k = 0:6
        f6sum_MX = f6sum_MX + ((alpha.MX.*r).^k)./factorial(k); 
        f6sum_MM = f6sum_MM + ((alpha.MM.*r).^k)./factorial(k);
        f6sum_XX = f6sum_XX + ((alpha.XX.*r).^k)./factorial(k);
    end
    f6_MX = f_r.*(1 - f6sum_MX.*exp(-alpha.MX.*r));
    f6_MM = f_r.*(1 - f6sum_MM.*exp(-alpha.MM.*r));
    f6_XX = f_r.*(1 - f6sum_XX.*exp(-alpha.XX.*r));

    % C8 damping functions
    f8sum_MX = 0;
    f8sum_MM = 0;
    f8sum_XX = 0;
    for k = 0:8
        f8sum_MX = f8sum_MX + ((alpha.MX.*r).^k)./factorial(k);
        f8sum_MM = f8sum_MM + ((alpha.MM.*r).^k)./factorial(k);
        f8sum_XX = f8sum_XX + ((alpha.XX.*r).^k)./factorial(k);
    end
    f8_MX = f_r.*(1 - f8sum_MX.*exp(-alpha.MX.*r));
    f8_MM = f_r.*(1 - f8sum_MM.*exp(-alpha.MM.*r));
    f8_XX = f_r.*(1 - f8sum_XX.*exp(-alpha.XX.*r));

    % Values of damping function at vdw cutoff
    f6_cutoff_MX = f_cutoff*(1 - sum(((alpha.MX.*RVDW_Cutoff).^(0:6))./factorial(0:6)).*exp(-alpha.MX.*RVDW_Cutoff));
    f6_cutoff_MM = f_cutoff*(1 - sum(((alpha.MM.*RVDW_Cutoff).^(0:6))./factorial(0:6)).*exp(-alpha.MM.*RVDW_Cutoff));
    f6_cutoff_XX = f_cutoff*(1 - sum(((alpha.XX.*RVDW_Cutoff).^(0:6))./factorial(0:6)).*exp(-alpha.XX.*RVDW_Cutoff));
    f8_cutoff_MX = f_cutoff*(1 - sum(((alpha.MX.*RVDW_Cutoff).^(0:8))./factorial(0:8)).*exp(-alpha.MX.*RVDW_Cutoff));
    f8_cutoff_MM = f_cutoff*(1 - sum(((alpha.MM.*RVDW_Cutoff).^(0:8))./factorial(0:8)).*exp(-alpha.MM.*RVDW_Cutoff));
    f8_cutoff_XX = f_cutoff*(1 - sum(((alpha.XX.*RVDW_Cutoff).^(0:8))./factorial(0:8)).*exp(-alpha.XX.*RVDW_Cutoff));

    % Calculate C6 damping derivatives
    df6sum_MX = 0;
    df6sum_MM = 0;
    df6sum_XX = 0;
    for k = 1:6
        df6sum_MX = df6sum_MX + k.*(alpha.MX.^k).*(r.^(k-1))./factorial(k);
        df6sum_MM = df6sum_MM + k.*(alpha.MM.^k).*(r.^(k-1))./factorial(k);
        df6sum_XX = df6sum_XX + k.*(alpha.XX.^k).*(r.^(k-1))./factorial(k);
    end
    df6_MX = f_r.*((alpha.MX.*exp(-alpha.MX.*r)).*f6sum_MX ...
        - (exp(-alpha.MX.*r)).*df6sum_MX) + df_r.*(1 - f6sum_MX.*exp(-alpha.MX.*r));
    df6_MM = f_r.*((alpha.MX.*exp(-alpha.MX.*r)).*f6sum_MM ...
        - (exp(-alpha.MM.*r)).*df6sum_MM) + df_r.*(1 - f6sum_MM.*exp(-alpha.MM.*r));
    df6_XX = f_r.*((alpha.MX.*exp(-alpha.MX.*r)).*f6sum_XX ...
        - (exp(-alpha.XX.*r)).*df6sum_XX) + df_r.*(1 - f6sum_XX.*exp(-alpha.XX.*r));

    %% Calculate C8 dispersion derivative with damping
    df8sum_MX = 0;
    df8sum_MM = 0;
    df8sum_XX = 0;
    for k = 1:8
        df8sum_MX = df8sum_MX + k.*(alpha.MX.^k).*(r.^(k-1))./factorial(k);
        df8sum_MM = df8sum_MM + k.*(alpha.MM.^k).*(r.^(k-1))./factorial(k);
        df8sum_XX = df8sum_XX + k.*(alpha.XX.^k).*(r.^(k-1))./factorial(k);
    end
    df8_MX = f_r.*((alpha.MX.*exp(-alpha.MX.*r)).*f8sum_MX ...
        - (exp(-alpha.MX.*r)).*df8sum_MX) + df_r.*(1 - f8sum_MX.*exp(-alpha.MX.*r));
    df8_MM = f_r.*((alpha.MX.*exp(-alpha.MX.*r)).*f8sum_MM ...
        - (exp(-alpha.MM.*r)).*df8sum_MM) + df_r.*(1 - f8sum_MM.*exp(-alpha.MM.*r));
    df8_XX = f_r.*((alpha.MX.*exp(-alpha.MX.*r)).*f8sum_XX ...
        - (exp(-alpha.XX.*r)).*df8sum_XX) + df_r.*(1 - f8sum_XX.*exp(-alpha.XX.*r));

%% Family of other damping functions related by a single formula, requiring van der waals radii
elseif C6Damping >= 3 && C6Damping <= 6
    
    % Define crystal radii source:
    % https://en.wikipedia.org/wiki/Ionic_radius
    R0.Li = 0.09; % nm
    R0.Na = 0.116; % nm
    R0.K  = 0.152; % nm
    R0.Rb = 0.166; % nm
    R0.Cs = 0.181; % nm
    R0.F  = 0.119; % nm
    R0.Cl = 0.167; % nm
    R0.Br = 0.182; % nm
    R0.I  = 0.206; % nm
    
    %% MMDRE Damping function. Citation:
    % "Transferable ab initio intermolecular potentials. 1. Derivation from methanol dimer and trimer calculations"
    % W.T.M. Mooij, F.B. van Duijneveldt, J.G.C.M. van Duijneveldt-van de Rijdt, B.P. van Eijck
    % J. Phys. Chem. A 1999, 103, 9872-9882.
    if C6Damping == 3
        a = -1;
        b = 7.19;
        m = 3;
        n = 2;
        
    %% PAMoC Damping function. Citation:
    % "Empirical correction to density functional theory for van der Waals interactions"
    % Q. Wu, W. Yang
    % J. Chem. Phys. 2002, 116, 515-524.
    elseif C6Damping == 4
        a = -1;
        b = 3.54;
        m = 3;
        n = 2;
        
    %% EHFSK Damping function. Citation:
    % "Hydrogen bonding and stacking interactions of nucleic acid base pairs: A density-functional-theory based treatment"
    % M. Elstner, P. Hobza, T. Frauenheim, S. Suhai, E. Kaxiras
    % J. Chem. Phys. 2001, 114, 5149-5155.
    elseif C6Damping == 5
        a = -1;
        b = 3;
        m = 7;
        n = 4;
    
    %% WY damping function. Citation:
    % "Empirical correction to density functional theory for van der Waals interactions"
    % Q. Wu, W. Yang
    % J. Chem. Phys. 2002, 116, 515-524.
    elseif C6Damping == 6
        a = exp(23);
        b = 23;
        m = 1;
        n = -1;
    end
    
    R0_MX = R0.(Metal) + R0.(Halide);
    R0_MM = R0.(Metal) + R0.(Metal);
    R0_XX = R0.(Halide) + R0.(Halide);
    
    % Damping functions
    f6_MX = f_r.*((1 + a.*exp(-b.*(r./R0_MX).^m) ).^n);
    f6_MM = f_r.*((1 + a.*exp(-b.*(r./R0_MM).^m) ).^n);
    f6_XX = f_r.*((1 + a.*exp(-b.*(r./R0_XX).^m) ).^n);
    
    f8_MX = f6_MX; % These are not separately defined
    f8_MM = f6_MM; % These are not separately defined
    f8_XX = f6_XX; % These are not separately defined

    % Derivative parts of damping functions
    df6_MX = f_r.*(-(a.*b.*m.*n.*(r./R0_MX).^(m - 1).*exp(-b.*(r./R0_MX).^m).*(a.*exp(-b.*(r/R0_MX).^m) + 1).^(n - 1))./R0_MX) ...
        + df_r.*((1 + a.*exp(-b.*(r./R0_MX).^m) ).^n);           
    df6_MM = f_r.*(-(a.*b.*m.*n.*(r./R0_MM).^(m - 1).*exp(-b.*(r./R0_MM).^m).*(a.*exp(-b.*(r/R0_MM).^m) + 1).^(n - 1))./R0_MM) ...
        + df_r.*((1 + a.*exp(-b.*(r./R0_MM).^m) ).^n);
    df6_XX = f_r.*(-(a.*b.*m.*n.*(r./R0_XX).^(m - 1).*exp(-b.*(r./R0_XX).^m).*(a.*exp(-b.*(r/R0_XX).^m) + 1).^(n - 1))./R0_XX) ...
        + df_r.*((1 + a.*exp(-b.*(r./R0_XX).^m) ).^n);
    
    df8_MX = df6_MX;
    df8_MM = df6_MM;
    df8_XX = df6_XX;
    
    % Values at vdw cutoff
    f6_cutoff_MX = f_cutoff.*((1 + a*exp(-b*(RVDW_Cutoff/R0_MX)^m) )^n);
    f6_cutoff_MM = f_cutoff.*((1 + a*exp(-b*(RVDW_Cutoff/R0_MM)^m) )^n);
    f6_cutoff_XX = f_cutoff.*((1 + a*exp(-b*(RVDW_Cutoff/R0_XX)^m) )^n);
    f8_cutoff_MX = f6_cutoff_MX;
    f8_cutoff_MM = f6_cutoff_MM;
    f8_cutoff_XX = f6_cutoff_XX;
end

%% Build PES: Plus-Minus (PM)

% Plus - Minus total potential
U_MetHal.Total = f_r.*k_0*(e_c^2)*q.(Metal)*q.(Halide)./(r) ...
    + B.MX*exp(-alpha.MX.*r) ...
    - f6_MX.*C.MX./(r.^6) ...
    - f8_MX.*D.MX./(r.^8);

% Components of Plus - Minus potential
U_MetHal.f = 1./r; % Electrostatics function f(r)
U_MetHal.g = - f6_MX.*C.MX./(r.^6) - f8_MX.*D.MX./(r.^8); % Dispersion g(r)
U_MetHal.h = B.MX*exp(-alpha.MX.*r) - k_0*(e_c^2)*q.(Metal)*q.(Halide)./(r) ...
    + f_r.*k_0*(e_c^2)*q.(Metal)*q.(Halide)./(r);% Short range repulsion h(r)

% Plus - Minus total derivative
U_MetHal.dTotal = -f_r.*k_0*(e_c^2)*q.(Metal)*q.(Halide)./(r.^2) ...
    + df_r.*k_0*(e_c^2)*q.(Metal)*q.(Halide)./(r) ...
    - alpha.MX*B.MX*exp(-alpha.MX.*r) ...
    + f6_MX.*6.*C.MX./(r.^7) ...
    - df6_MX.*C.MX./(r.^6) ...
    + f8_MX.*8.*D.MX./(r.^9) ...
    - df8_MX.*D.MX./(r.^8);

% Components of Plus - Minus forces
U_MetHal.df = 1./(r.^2);% Electrostatics function -df(r)/dr
U_MetHal.dg = - f6_MX.*6.*C.MX./(r.^7) + df6_MX.*C.MX./(r.^6) ...
    - f8_MX.*8.*D.MX./(r.^9) + df8_MX.*D.MX./(r.^8); % Dispersion -dg(r)/dr
U_MetHal.dh = alpha.MX*B.MX*exp(-alpha.MX.*r) - k_0*(e_c^2)*q.(Metal)*q.(Halide)./(r.^2) ...
	+ f_r.*k_0*(e_c^2)*q.(Metal)*q.(Halide)./(r.^2) ...
    - df_r.*k_0*(e_c^2)*q.(Metal)*q.(Halide)./(r);% Short range repulsion -dh(r)/dr

% Shift the Plus - Minus potential
if contains(vdw_modifier,'potential-shift','IgnoreCase',true)
    EVDW_Cutoff = B.MX*exp(-alpha.MX*RVDW_Cutoff) ...
        - k_0*(e_c^2)*q.(Metal)*q.(Halide)/(RVDW_Cutoff) ...
        + f_cutoff*k_0*(e_c^2)*q.(Metal)*q.(Halide)/(RVDW_Cutoff) ...
        - f6_cutoff_MX*C.MX/(RVDW_Cutoff^6) ...
        - f8_cutoff_MX*D.MX/(RVDW_Cutoff^8);
    
    % Shift by the dispersion energy at vdw cutoff radius. only affects one
    % energy component, not derivatives (i.e. forces)
    U_MetHal.Total = U_MetHal.Total - EVDW_Cutoff;
    U_MetHal.g = U_MetHal.g - EVDW_Cutoff;
end

% Remove infinities (replace with zero)
U_MetHal = Remove_Infinities(U_MetHal);

%% Build PES: Plus - Plus (PP)

% Plus - Plus Total Potential
U_MetMet.Total = k_0*(e_c^2)*q.(Metal).*q.(Metal)./(r) ...
    + B.MM*exp(-alpha.MM.*r) ...
    - f6_MM.*C.MM./(r.^6) ...
    - f8_MM.*D.MM./(r.^8);

% Components of Plus - Plus potential
U_MetMet.f = 1./r; % Electrostatics function f(r)
U_MetMet.g = - f6_MM.*C.MM./(r.^6) - f8_MM.*D.MM./(r.^8); % Dispersion g(r)
U_MetMet.h = B.MM*exp(-alpha.MM.*r); % Short range repulsion

% Plus - Plus Total Derivative
U_MetMet.dTotal = -k_0*(e_c^2)*q.(Metal).*q.(Metal)./(r.^2) ...
    - alpha.MM*B.MM*exp(-alpha.MM.*r)...
    + f6_MM.*6.*C.MM./(r.^7) ...
    - df6_MM.*C.MM./(r.^6) ...
    + f8_MM.*8.*D.MM./(r.^9) ...
    - df8_MM.*D.MM./(r.^8);

% Components of Plus - Plus forces
U_MetMet.df = 1./(r.^2); % Electrostatics function f(r)
U_MetMet.dg = - f6_MM.*6.*C.MM./(r.^7) + df6_MM.*C.MM./(r.^6) ...
    - f8_MM.*8.*D.MM./(r.^9) + df8_MM.*D.MM./(r.^8); % Dispersion g(r)
U_MetMet.dh = alpha.MM*B.MM*exp(-alpha.MM.*r); % Short range repulsion

% Shift the Plus - Plus potential
if contains(vdw_modifier,'potential-shift','IgnoreCase',true)
    EVDW_Cutoff = B.MM*exp(-alpha.MM*RVDW_Cutoff) ...
    - f6_cutoff_MM*C.MM/(RVDW_Cutoff^6) ...
    - f8_cutoff_MM*D.MM/(RVDW_Cutoff^8);
    
    % Shift by the dispersion energy at vdw cutoff radius. only affects one
    % energy component, not derivatives (i.e. forces)
    U_MetMet.Total = U_MetMet.Total - EVDW_Cutoff;
    U_MetMet.g = U_MetMet.g - EVDW_Cutoff;
end

% Remove infinities (replace with zero)
U_MetMet = Remove_Infinities(U_MetMet);

%% Build Model Table: Minus - Minus (MM)
% Minus - Minus Total Potential
U_HalHal.Total = k_0*(e_c^2)*q.(Halide).*q.(Halide)./(r)...
    + B.XX*exp(-alpha.XX.*r) ...
    - f6_XX.*C.XX./(r.^6)...
    - f8_XX.*D.XX./(r.^8);

% Components of Minus - Minus potential
U_HalHal.f = 1./r; % Electrostatics function f(r)
U_HalHal.g = - f6_XX.*C.XX./(r.^6) - f8_XX.*D.XX./(r.^8); % Dispersion g(r)
U_HalHal.h = B.XX*exp(-alpha.XX.*r); % Short range repulsion

% Minus - Minus Total Derivative
U_HalHal.dTotal = -k_0*(e_c^2)*q.(Halide).*q.(Halide)./(r.^2)...
    - alpha.XX*B.XX*exp(-alpha.XX.*r)...
    + f6_XX.*6.*C.XX./(r.^7)...
    - df6_XX.*C.XX./(r.^6) ...
    + f8_XX.*8.*D.XX./(r.^9) ...
    - df8_XX.*D.XX./(r.^8);

% Components of Minus - Minus forces
U_HalHal.df = 1./(r.^2); % Electrostatics function f(r)
U_HalHal.dg = - f6_XX.*6.*C.XX./(r.^7) + df6_XX.*C.XX./(r.^6) ...
    - f8_XX.*8.*D.XX./(r.^9) + df8_XX.*D.XX./(r.^8); % Dispersion g(r)
U_HalHal.dh = alpha.XX*B.XX*exp(-alpha.XX.*r); % Short range repulsion

% Shift the Minus - Minus potential
if contains(vdw_modifier,'potential-shift','IgnoreCase',true)
    EVDW_Cutoff = B.XX*exp(-alpha.XX*RVDW_Cutoff) ...
    - f6_cutoff_XX*C.XX/(RVDW_Cutoff^6)...
    - f8_cutoff_XX*D.XX/(RVDW_Cutoff^8);
    
    % Shift by the dispersion energy at vdw cutoff radius. only affects one
    % energy component, not derivatives (i.e. forces)
    U_HalHal.Total = U_HalHal.Total - EVDW_Cutoff;
    U_HalHal.g = U_HalHal.g - EVDW_Cutoff;
end

% remove infinities
U_HalHal = Remove_Infinities(U_HalHal);

%% Print Potentials
PM = [r ; U_MetHal.f ; U_MetHal.df ; U_MetHal.g ; U_MetHal.dg ; U_MetHal.h ; U_MetHal.dh];
U_PM_out = deblank( sprintf(['%16.10e   %16.10e %16.10e   %16.10e %16.10e   %16.10e %16.10e' newline],PM(:)) );

PP = [r ; U_MetMet.f ; U_MetMet.df ; U_MetMet.g ; U_MetMet.dg ; U_MetMet.h ; U_MetMet.dh];
U_PP_out = deblank( sprintf(['%16.10e   %16.10e %16.10e   %16.10e %16.10e   %16.10e %16.10e' newline],PP(:)) );

MM = [r ; U_HalHal.f ; U_HalHal.df ; U_HalHal.g ; U_HalHal.dg ; U_HalHal.h ; U_HalHal.dh];
U_MM_out = deblank( sprintf(['%16.10e   %16.10e %16.10e   %16.10e %16.10e   %16.10e %16.10e' newline],MM(:)) );

%% Plot potentials if plotswitch chosen
if plotswitch
    figure;
    
    % Options
    lw=2;
    fs=25;

    h = cell(1,8);
    hold on

    h{1} = plot(r,-U_MetHal.dg,'Color','b','LineWidth',lw,'LineStyle','-');
    h{2} = plot(r,-U_MetMet.dg,'Color','r','LineWidth',lw,'LineStyle','-');
    h{3} = plot(r,-U_HalHal.dg,'Color','g','LineWidth',lw,'LineStyle','-');
    
    dr = r(2) - r(1);
    r_dr = r(2:end) - dr/2;
    
    h{1} = plot(r_dr,diff(U_MetHal.g)./dr,'Color','k','LineWidth',lw,'LineStyle',':');
    h{2} = plot(r_dr,diff(U_MetMet.g)./dr,'Color','k','LineWidth',lw,'LineStyle',':');
    h{3} = plot(r_dr,diff(U_HalHal.g)./dr,'Color','k','LineWidth',lw,'LineStyle',':');

    title(['Plot of TF Potentials for ' Salt],...
        'Interpreter','latex','fontsize',fs)

    set(gca,'box','on','TickLabelInterpreter','latex');
    set(gca,'XMinorTick','on','YMinorTick','on','FontSize',fs);
    xlabel('r (nm)','fontsize',fs,'Interpreter','latex');
    ylabel('Pair Potential (kJ mol$^{-1}$)','fontsize',fs,'Interpreter','latex');

    ylim([-1000 1000]);
    xlim([Startpoint Endpoint]);

    % Blank line
    hline = refline([0 0]);
    hline.Color = 'k';
    hline.LineWidth = lw-1;
    hline.LineStyle = '--';
    leg1 = legend([h{:}],{['TF - ' Salt] ['TF - ' Metal Metal] ['TF - ' Halide Halide]});
    leg1.Interpreter = 'latex';
end

end