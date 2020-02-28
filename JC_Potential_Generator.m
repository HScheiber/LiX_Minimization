% Generates Jung-Cheatham pair potential energy tables from 'Startpoint' up to a given
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

function [U_PM_out, U_PP_out, U_MM_out] = JC_Potential_Generator(Startpoint,...
    Endpoint,Spacing,Salt,Parameters,plotswitch,vdw_modifier,RVDW_Cutoff,CRDamping,C6Damping)

%% Conversion factors and fundamental constants
nm_per_m = 1e+9; % nm per m
NA = 6.0221409e23; % Molecules per mole
e_c = 1.60217662e-19; % Elementary charge in Coulombs
epsilon_0 = (8.854187817620e-12)*1000/(nm_per_m*NA); % Vacuum Permittivity C^2 mol kJ^-1 nm^-1
k_0 = 1/(4*pi*epsilon_0); % Coulomb constant in kJ nm C^-2 mol^-1
nm_per_Ang = 0.1; % nm per Angstrom

%% Close-Range Sigmoid Damping Parameters
r_d = 0.10; % If close-range dispersion damping is on, this is the value of the sigmoid's midpoint
b = 150; % sigmoid "steepness" for damping

[Metal,Halide] = Separate_Metal_Halide(Salt);

%% Parameter: q (charge)
q.Li =  1; % atomic
q.Na =  1; % atomic
q.K  =  1; % atomic
q.Rb =  1; % atomic
q.Cs =  1; % atomic

q.F  = -1; % atomic
q.Cl = -1; % atomic
q.Br = -1; % atomic
q.I  = -1; % atomic


%% Parameters from input for LJ potential
sigma_MM = Parameters(1,1);
sigma_XX = Parameters(1,2);
sigma_MX = Parameters(1,3);

epsilon_MM = Parameters(2,1);
epsilon_XX = Parameters(2,2);
epsilon_MX = Parameters(2,3);

% Change parameters into A/r12 - B/r6 format
A_MM = 4*epsilon_MM*(sigma_MM^12);
B_MM = 4*epsilon_MM*(sigma_MM^6);

A_XX = 4*epsilon_XX*(sigma_XX^12);
B_XX = 4*epsilon_XX*(sigma_XX^6);

A_MX = 4*epsilon_MX*(sigma_MX^12);
B_MX = 4*epsilon_MX*(sigma_MX^6);

%% Generate range (r) in nm
r = Startpoint:Spacing:Endpoint;

%% If Damping at close range
if CRDamping
    f_r = 1./(1 + exp(-b.*(r - r_d))); % sigmoid damping function
    df_r = (b.*exp(-b.*(r - r_d)))./((1 + exp(-b.*(r - r_d))).^2); % sigmoid damping function derivative
    f_cutoff = 1./(1 + exp(-b.*(RVDW_Cutoff - r_d)));  % value of f at vdw cutoff
else
	f_r = 1; % No damping
    df_r = 0; % Damping derivative is zero
    f_cutoff = 1; % value of f at vdw cutoff
end

%% Medium range dispersion damping: No Damping (default)
if C6Damping == 0
    % No-damping function
    f6_MX = f_r;
    f6_MM = f_r;
    f6_XX = f_r;

    % Derivative parts of no-damping function   
    df6_MX = df_r;
    df6_MM = df_r;
    df6_XX = df_r;

    % Values at vdw cutoff
    f6_cutoff_MX = f_cutoff;
    f6_cutoff_MM = f_cutoff;
    f6_cutoff_XX = f_cutoff;
    
%% BJ-Type Rational Damping
elseif C6Damping == 1
    
    % Generate conversion factors
    Bohr_nm = 0.0529177; % a_0 - > Angstrom
    c6conv = 1e-3/2625.4999/((0.052917726)^6); % J/mol nm^6 - > au (from sourcecode)
    J_kJ = 1e-3; % J - > kJ
    Ha_kJmol = 2625.4999; % Ha - > kJ/mol
    c6units = (1/c6conv)*J_kJ; % au - > kJ/mol nm^6
    c8units = (Ha_kJmol)*(Bohr_nm^8); % au - > kJ/mol nm^8
    
    % Factor used to calculate C8
    sqrt_Q.Li = 5.019869340000000;
    sqrt_Q.Na = 6.585855360000000;
    sqrt_Q.K  = 7.977627530000000;
    sqrt_Q.Rb = 9.554616980000000;
    sqrt_Q.Cs = 11.02204549000000;
    sqrt_Q.F  = 2.388252500000000;
    sqrt_Q.Cl = 3.729323560000000;
    sqrt_Q.Br = 4.590896470000000;
    sqrt_Q.I  = 5.533218150000000;
    
    % Calculate C8 (needed for cutoff radius)
    C8_MX = 3.0*(B_MX/c6units)*sqrt_Q.(Metal)*sqrt_Q.(Halide)*c8units; % in kJ/mol nm^8
    C8_MM = 3.0*(B_MM/c6units)*sqrt_Q.(Metal)*sqrt_Q.(Metal)*c8units; % in kJ/mol nm^8
    C8_XX = 3.0*(B_XX/c6units)*sqrt_Q.(Halide)*sqrt_Q.(Halide)*c8units; % in kJ/mol nm^8
    
    % Damping distances (no C8 term so define wrt crystal radii)
    R0_MX = sqrt(C8_MX/B_MX); % in Angstroms
    R0_MM = sqrt(C8_MM/B_MM); % in Angstroms  
    R0_XX = sqrt(C8_XX/B_XX); % in Angstroms
    
    % Damping functions (unitless)
    f6_MX = f_r./( 1 + ( R0_MX ./ r ).^6 );
    f6_MM = f_r./( 1 + ( R0_MM ./ r ).^6 );
    f6_XX = f_r./( 1 + ( R0_XX ./ r ).^6 );
    
    % Values of damping function at vdw cutoff
    f6_cutoff_MX = f_cutoff/( 1 + ( R0_MX / RVDW_Cutoff )^6 );
    f6_cutoff_MM = f_cutoff/( 1 + ( R0_MM / RVDW_Cutoff )^6 );
    f6_cutoff_XX = f_cutoff/( 1 + ( R0_XX / RVDW_Cutoff )^6 );
    
    % Derivative of damping functions
    df6_MX = f_r.*(6.*(R0_MX.^6).*(r.^5)./(((r.^6) + (R0_MX.^6)).^2)) + df_r./( 1 + ( R0_MX ./ r ).^6 );
    df6_MM = f_r.*(6.*(R0_MM.^6).*(r.^5)./(((r.^6) + (R0_MM.^6)).^2)) + df_r./( 1 + ( R0_MM ./ r ).^6 );
    df6_XX = f_r.*(6.*(R0_XX.^6).*(r.^5)./(((r.^6) + (R0_XX.^6)).^2)) + df_r./( 1 + ( R0_XX ./ r ).^6 );
    
%% Tang and Toennies Damping function. Cite:
% "An improved simple model for the van der Waals potential based on universal damping functions for the dispersion coefficients."
% K. T. Tang, J. P. Toennies
% J. Chem. Phys. 1984, 80, 3726-3741.
elseif C6Damping == 2

    % Use the TF hardness parameters (no equivalent parameter for JC!)
    rho.LiF = 0.299*nm_per_Ang; % nm
    rho.LiCl = 0.342*nm_per_Ang; % nm
    rho.LiBr = 0.353*nm_per_Ang; % nm
    rho.LiI = 0.430*nm_per_Ang; % nm

    rho.NaF = 0.330*nm_per_Ang; % nm
    rho.NaCl = 0.317*nm_per_Ang; % nm
    rho.NaBr = 0.340*nm_per_Ang; % nm
    rho.NaI = 0.386*nm_per_Ang; % nm

    rho.KF = 0.338*nm_per_Ang; % nm
    rho.KCl = 0.337*nm_per_Ang; % nm
    rho.KBr = 0.335*nm_per_Ang; % nm
    rho.KI = 0.355*nm_per_Ang; % nm

    rho.RbF = 0.328*nm_per_Ang; % nm
    rho.RbCl = 0.318*nm_per_Ang; % nm
    rho.RbBr = 0.335*nm_per_Ang; % nm
    rho.RbI = 0.337*nm_per_Ang; % nm

    rho.CsF = 0.282*nm_per_Ang; % nm
    alpha = rho.(Salt);
    
    % C6 damping functions
    f6sum_MX = 0;
    f6sum_MM = 0;
    f6sum_XX = 0;
    for k = 0:6
        f6sum_MX = f6sum_MX + ((alpha.*r).^k)./factorial(k); 
        f6sum_MM = f6sum_MM + ((alpha.*r).^k)./factorial(k);
        f6sum_XX = f6sum_XX + ((alpha.*r).^k)./factorial(k);
    end
    f6_MX = f_r.*(1 - f6sum_MX.*exp(-alpha.*r));
    f6_MM = f_r.*(1 - f6sum_MM.*exp(-alpha.*r));
    f6_XX = f_r.*(1 - f6sum_XX.*exp(-alpha.*r));

    % Calculate C6 damping derivatives
    df6sum_MX = 0;
    df6sum_MM = 0;
    df6sum_XX = 0;
    for k = 1:6
        df6sum_MX = df6sum_MX + k.*(alpha.^k).*(r.^(k-1))./factorial(k);
        df6sum_MM = df6sum_MM + k.*(alpha.^k).*(r.^(k-1))./factorial(k);
        df6sum_XX = df6sum_XX + k.*(alpha.^k).*(r.^(k-1))./factorial(k);
    end
    df6_MX = f_r.*((alpha.*exp(-alpha.*r)).*f6sum_MX ...
        - (exp(-alpha.*r)).*df6sum_MX) + df_r.*(1 - f6sum_MX.*exp(-alpha.*r));
    df6_MM = f_r.*((alpha.*exp(-alpha.*r)).*f6sum_MM ...
        - (exp(-alpha.*r)).*df6sum_MM) + df_r.*(1 - f6sum_MM.*exp(-alpha.*r));
    df6_XX = f_r.*((alpha.*exp(-alpha.*r)).*f6sum_XX ...
        - (exp(-alpha.*r)).*df6sum_XX) + df_r.*(1 - f6sum_XX.*exp(-alpha.*r));
    
    % Values of damping function at vdw cutoff
    f6_cutoff_MX = f_cutoff*(1 - sum(((alpha.*RVDW_Cutoff).^(0:6))./factorial(0:6)).*exp(-alpha.*RVDW_Cutoff));
    f6_cutoff_MM = f_cutoff*(1 - sum(((alpha.*RVDW_Cutoff).^(0:6))./factorial(0:6)).*exp(-alpha.*RVDW_Cutoff));
    f6_cutoff_XX = f_cutoff*(1 - sum(((alpha.*RVDW_Cutoff).^(0:6))./factorial(0:6)).*exp(-alpha.*RVDW_Cutoff));

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

    % Derivative parts of damping functions
    df6_MX = f_r.*(-(a.*b.*m.*n.*(r./R0_MX).^(m - 1).*exp(-b.*(r./R0_MX).^m).*(a.*exp(-b.*(r/R0_MX).^m) + 1).^(n - 1))./R0_MX) ...
        + df_r.*((1 + a.*exp(-b.*(r./R0_MX).^m) ).^n);            
    df6_MM = f_r.*(-(a.*b.*m.*n.*(r./R0_MM).^(m - 1).*exp(-b.*(r./R0_MM).^m).*(a.*exp(-b.*(r/R0_MM).^m) + 1).^(n - 1))./R0_MM) ...
        + df_r.*((1 + a.*exp(-b.*(r./R0_MM).^m) ).^n);
    df6_XX = f_r.*(-(a.*b.*m.*n.*(r./R0_XX).^(m - 1).*exp(-b.*(r./R0_XX).^m).*(a.*exp(-b.*(r/R0_XX).^m) + 1).^(n - 1))./R0_XX) ...
        + df_r.*((1 + a.*exp(-b.*(r./R0_XX).^m) ).^n);
    
    % Values at vdw cutoff
    f6_cutoff_MX = f_cutoff*((1 + a*exp(-b*(RVDW_Cutoff/R0_MX)^m) )^n);
    f6_cutoff_MM = f_cutoff*((1 + a*exp(-b*(RVDW_Cutoff/R0_MM)^m) )^n);
    f6_cutoff_XX = f_cutoff*((1 + a*exp(-b*(RVDW_Cutoff/R0_XX)^m) )^n);
end

%% Build PES: Plus-Minus
% Plus - Minus total potential
U_MetHal.Total = f_r.*k_0*(e_c^2)*q.(Metal)*q.(Halide)./r ...
+ A_MX./(r.^12) ...
- f6_MX.*B_MX./(r.^6);

% Components of Plus - Minus potential
U_MetHal.f = 1./r; % Electrostatics function f(r)
U_MetHal.g = -f6_MX.*B_MX./(r.^6); % Dispersion g(r)
U_MetHal.h = A_MX./(r.^12) - k_0*(e_c^2)*q.(Metal)*q.(Halide)./(r) + ...
    f_r.*k_0*(e_c^2)*q.(Metal)*q.(Halide)./(r); % Short range repulsion

% Plus - Minus total derivative
U_MetHal.dTotal = -f_r.*k_0*(e_c^2)*q.(Metal)*q.(Halide)./(r.^2) ...
    + df_r.*k_0*(e_c^2)*q.(Metal)*q.(Halide)./(r)...
    - A_MX.*12./(r.^13) ...
    + f6_MX.*B_MX.*6./(r.^7) ...
    - df6_MX.*B_MX./(r.^6);

% Components of Plus - Minus forces
U_MetHal.df = 1./(r.^2);% Electrostatics function
U_MetHal.dg = - f6_MX.*B_MX.*6./(r.^7) + df6_MX.*B_MX./(r.^6); % Dispersion
U_MetHal.dh = + A_MX.*12./(r.^13) - k_0*(e_c^2)*q.(Metal)*q.(Halide)./(r.^2) ...
    + f_r.*k_0*(e_c^2)*q.(Metal)*q.(Halide)./(r.^2) ...
    - df_r.*k_0*(e_c^2)*q.(Metal)*q.(Halide)./(r);% Short range repulsion

% Shift the Plus - Minus potential
if contains(vdw_modifier,'potential-shift','IgnoreCase',true)
    EVDW_Cutoff = A_MX/(RVDW_Cutoff^12) ...
        - k_0*(e_c^2)*q.(Metal)*q.(Halide)/(RVDW_Cutoff) ...
        + f_cutoff*k_0*(e_c^2)*q.(Metal)*q.(Halide)/(RVDW_Cutoff) ...
        - f6_cutoff_MX*B_MX/(RVDW_Cutoff^6);
    
    % Shift by the dispersion energy at vdw cutoff radius. only affects one
    % energy component, not derivatives (i.e. forces)
    U_MetHal.Total = U_MetHal.Total - EVDW_Cutoff;
    U_MetHal.g = U_MetHal.g - EVDW_Cutoff;
end

% remove infinities
U_MetHal = Remove_Infinities(U_MetHal);

%% Build PES: Plus - Plus
% Plus - Plus total potential
U_MetMet.Total = k_0*(e_c^2)*q.(Metal)*q.(Metal)./(r) ...
    + A_MM./(r.^12) ...
    - f6_MM.*B_MM./(r.^6);

% Plus - Plus potential components
U_MetMet.f = 1./r; % Electrostatics function f(r)
U_MetMet.g = -f6_MM.*B_MM./(r.^6); % Dispersion g(r)
U_MetMet.h = A_MM./(r.^12);% Short range repulsion

% Plus - Plus total derivative
U_MetMet.dTotal = -k_0*(e_c^2)*q.(Metal)*q.(Metal)./(r.^2) ...
    - A_MM.*12./(r.^13) ...
    + f6_MM.*B_MM.*6./(r.^7) ...
    - df6_MM.*B_MM./(r.^6);

% Components of Plus - Plus forces 
U_MetMet.df = 1./(r.^2);% Electrostatics function
U_MetMet.dg = - f6_MM.*B_MM.*6./(r.^7) + df6_MM.*B_MM./(r.^6); % Dispersion
U_MetMet.dh = + A_MM.*12./(r.^13);% Short range repulsion

% Shift the Plus - Plus potential
if contains(vdw_modifier,'potential-shift','IgnoreCase',true)
    EVDW_Cutoff = A_MM./(RVDW_Cutoff.^12) ...
    - f6_cutoff_MM*B_MM./(RVDW_Cutoff.^6);
    
    % Shift by the dispersion energy at vdw cutoff radius. only affects one
    % energy component, not derivatives (i.e. forces)
    U_MetMet.Total = U_MetMet.Total - EVDW_Cutoff;
    U_MetMet.g = U_MetMet.g - EVDW_Cutoff;
end

% remove infinities
U_MetMet = Remove_Infinities(U_MetMet);

%% Build PES: Minus - Minus
% Minus - Minus total potential
U_HalHal.Total = k_0*(e_c^2)*q.(Halide)*q.(Halide)./(r) ...
    + A_XX./(r.^12) ...
    - f6_XX.*B_XX./(r.^6);

% Minus - Minus potential components
U_HalHal.f = 1./r; % Electrostatics function f(r)
U_HalHal.g = -f6_XX.*B_XX./(r.^6); % Dispersion g(r)
U_HalHal.h = A_XX./(r.^12);% Short range repulsion

% Minus - Minus total derivative
U_HalHal.dTotal = -k_0*(e_c^2)*q.(Halide)*q.(Halide)./(r.^2) ...
    - A_XX.*12./(r.^13) ...
    + f6_XX.*B_XX.*6./(r.^7) ...
    - df6_XX.*B_XX./(r.^6);

% Components of Minus - Minus forces
U_HalHal.df = 1./(r.^2);% Electrostatics function
U_HalHal.dg = - f6_XX.*B_XX.*6./(r.^7) + df6_XX.*B_XX./(r.^6); % Dispersion
U_HalHal.dh = + A_XX.*12./(r.^13);% Short range repulsion

% Shift Minus - Minus potential
if contains(vdw_modifier,'potential-shift','IgnoreCase',true)
    EVDW_Cutoff = A_XX./(RVDW_Cutoff.^12) ...
    - f6_cutoff_XX*B_XX./(RVDW_Cutoff.^6);
    
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
if plotswitch % For testing
    figure;
    % Options
    lw=2;
    fs=25;

    h = cell(1,8);
    hold on

    h{1} = plot(r,U_MetHal.g,'Color','b','LineWidth',lw,'LineStyle','-');
    h{2} = plot(r,U_MetMet.g,'Color','r','LineWidth',lw,'LineStyle','-');
    h{3} = plot(r,U_HalHal.g,'Color','g','LineWidth',lw,'LineStyle','-');

%     dr = r(2) - r(1);
%     r_dr = r(2:end) - dr/2;
%     
%     h{4} = plot(r_dr,diff(U_MetHal.g)./dr,'Color','k','LineWidth',lw,'LineStyle',':');
%     h{5} = plot(r_dr,diff(U_MetMet.g)./dr,'Color','k','LineWidth',lw,'LineStyle',':');
%     h{6} = plot(r_dr,diff(U_HalHal.g)./dr,'Color','k','LineWidth',lw,'LineStyle',':');
    
    title(['Plot of JC Potentials for ' Salt],...
        'Interpreter','latex','fontsize',fs)

    set(gca,'box','on','TickLabelInterpreter','latex');
    set(gca,'XMinorTick','on','YMinorTick','on','FontSize',fs);
    xlabel('r (nm)','fontsize',fs,'Interpreter','latex');
    ylabel('Pair Potential (kJ mol$^{-1}$)','fontsize',fs,'Interpreter','latex');

    ylim([-100 50]);
    xlim([Startpoint Endpoint]);

    % Blank line
    hline = refline([0 0]);
    hline.Color = 'k';
    hline.LineWidth = lw-1;
    hline.LineStyle = '--';
    leg1 = legend([h{:}],{['JC - ' Salt] ['JC - ' Metal Metal] ['JC - ' Halide Halide]});
    leg1.Interpreter = 'latex';
end

end