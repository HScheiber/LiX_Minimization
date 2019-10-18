% Generates pair potential energy surfaces from 'Startpoint' up to a given
% length input 'Endpoint'. Plotswitch is a logical true or false to
% determine whether to plot the resulting potentials
% Recommended spacing is 0.005 angstroms or 0.0005 nm
% INPUT UNITS MUST BE ALL IN NANOMETERS. OUTPUTS ARE IN NANOMETERS AND
% kJ/mol
% DS is the dispersion scaling factor
function [U_PM_out, U_PP_out, U_MM_out] = TF_Potential_Generator(Startpoint,Endpoint,Spacing,Salt,plotswitch,DS)
%% Conversion factors and fundamental constants
kj_per_erg = 1e-10; % kJ per erg
Ang_per_cm = 1e+8; % Angstroms per cm
nm_per_cm = 1e+7; % nm per cm
Ang_per_m = 1e+10; % Angstroms per m
nm_per_m = 1e+9; % nm per m
NA = 6.0221409e23; % Molecules per mole
e_c = 1.60217662e-19; % Elementary charge in Coulombs
epsilon_0 = (8.854187817620e-12)*1000/(nm_per_m*NA); % Vacuum Permittivity C^2 mol kJ^-1 nm^-1
k_0 = 1/(4*pi*epsilon_0); % Coulomb constant in kJ nm C^-2 mol^-1
C_unit = 1e-60; % erg cm^6
D_unit = 1e-76; % erg cm^8
nm_per_Ang = 0.1; % nm per Angstrom

%% Huggins-Mayer Dipole-Dipole Dispersion Parameter C: PM = +-   PP = ++     MM = --
C.LiF.PP  = (0.073*C_unit)*kj_per_erg*(nm_per_cm^6)*NA; % kJ/mol nm^6
C.LiF.PM  = (0.8*C_unit)*kj_per_erg*(nm_per_cm^6)*NA; % kJ/mol nm^6
C.LiF.MM  = (14.5*C_unit)*kj_per_erg*(nm_per_cm^6)*NA; % kJ/mol nm^6

C.LiCl.PP = (0.073*C_unit)*kj_per_erg*(nm_per_cm^6)*NA; % kJ/mol nm^6
C.LiCl.PM = (2.0*C_unit)*kj_per_erg*(nm_per_cm^6)*NA; % kJ/mol nm^6
C.LiCl.MM = (111*C_unit)*kj_per_erg*(nm_per_cm^6)*NA; % kJ/mol nm^6

C.LiBr.PP = (0.073*C_unit)*kj_per_erg*(nm_per_cm^6)*NA; % kJ/mol nm^6
C.LiBr.PM = (2.5*C_unit)*kj_per_erg*(nm_per_cm^6)*NA; % kJ/mol nm^6
C.LiBr.MM = (185*C_unit)*kj_per_erg*(nm_per_cm^6)*NA; % kJ/mol nm^6

C.LiI.PP = (0.073*C_unit)*kj_per_erg*(nm_per_cm^6)*NA; % kJ/mol nm^6
C.LiI.PM = (3.3*C_unit)*kj_per_erg*(nm_per_cm^6)*NA; % kJ/mol nm^6
C.LiI.MM = (378*C_unit)*kj_per_erg*(nm_per_cm^6)*NA; % kJ/mol nm^6

C.NaF.PP  = (1.68*C_unit)*kj_per_erg*(nm_per_cm^6)*NA; % kJ/mol nm^6
C.NaF.PM  = (4.5*C_unit)*kj_per_erg*(nm_per_cm^6)*NA; % kJ/mol nm^6
C.NaF.MM  = (16.5*C_unit)*kj_per_erg*(nm_per_cm^6)*NA; % kJ/mol nm^6

C.NaCl.PP  = (1.68*C_unit)*kj_per_erg*(nm_per_cm^6)*NA; % kJ/mol nm^6
C.NaCl.PM  = (11.2*C_unit)*kj_per_erg*(nm_per_cm^6)*NA; % kJ/mol nm^6
C.NaCl.MM  = (116*C_unit)*kj_per_erg*(nm_per_cm^6)*NA; % kJ/mol nm^6

C.NaBr.PP  = (1.68*C_unit)*kj_per_erg*(nm_per_cm^6)*NA; % kJ/mol nm^6
C.NaBr.PM  = (14.0*C_unit)*kj_per_erg*(nm_per_cm^6)*NA; % kJ/mol nm^6
C.NaBr.MM  = (196*C_unit)*kj_per_erg*(nm_per_cm^6)*NA; % kJ/mol nm^6

C.NaI.PP  = (1.68*C_unit)*kj_per_erg*(nm_per_cm^6)*NA; % kJ/mol nm^6
C.NaI.PM  = (19.1*C_unit)*kj_per_erg*(nm_per_cm^6)*NA; % kJ/mol nm^6
C.NaI.MM  = (392*C_unit)*kj_per_erg*(nm_per_cm^6)*NA; % kJ/mol nm^6

C.KF.PP  = (24.3*C_unit)*kj_per_erg*(nm_per_cm^6)*NA; % kJ/mol nm^6
C.KF.PM  = (19.5*C_unit)*kj_per_erg*(nm_per_cm^6)*NA; % kJ/mol nm^6
C.KF.MM  = (18.6*C_unit)*kj_per_erg*(nm_per_cm^6)*NA; % kJ/mol nm^6

C.KCl.PP  = (24.3*C_unit)*kj_per_erg*(nm_per_cm^6)*NA; % kJ/mol nm^6
C.KCl.PM  = (48*C_unit)*kj_per_erg*(nm_per_cm^6)*NA; % kJ/mol nm^6
C.KCl.MM  = (124.5*C_unit)*kj_per_erg*(nm_per_cm^6)*NA; % kJ/mol nm^6

C.KBr.PP  = (24.3*C_unit)*kj_per_erg*(nm_per_cm^6)*NA; % kJ/mol nm^6
C.KBr.PM  = (60*C_unit)*kj_per_erg*(nm_per_cm^6)*NA; % kJ/mol nm^6
C.KBr.MM  = (206*C_unit)*kj_per_erg*(nm_per_cm^6)*NA; % kJ/mol nm^6

C.KI.PP  = (24.3*C_unit)*kj_per_erg*(nm_per_cm^6)*NA; % kJ/mol nm^6
C.KI.PM  = (82*C_unit)*kj_per_erg*(nm_per_cm^6)*NA; % kJ/mol nm^6
C.KI.MM  = (403*C_unit)*kj_per_erg*(nm_per_cm^6)*NA; % kJ/mol nm^6

C.RbF.PP  = (59.4*C_unit)*kj_per_erg*(nm_per_cm^6)*NA; % kJ/mol nm^6
C.RbF.PM  = (31*C_unit)*kj_per_erg*(nm_per_cm^6)*NA; % kJ/mol nm^6
C.RbF.MM  = (18.9*C_unit)*kj_per_erg*(nm_per_cm^6)*NA; % kJ/mol nm^6

C.RbCl.PP  = (59.4*C_unit)*kj_per_erg*(nm_per_cm^6)*NA; % kJ/mol nm^6
C.RbCl.PM  = (79*C_unit)*kj_per_erg*(nm_per_cm^6)*NA; % kJ/mol nm^6
C.RbCl.MM  = (130*C_unit)*kj_per_erg*(nm_per_cm^6)*NA; % kJ/mol nm^6

C.RbBr.PP  = (59.4*C_unit)*kj_per_erg*(nm_per_cm^6)*NA; % kJ/mol nm^6
C.RbBr.PM  = (99*C_unit)*kj_per_erg*(nm_per_cm^6)*NA; % kJ/mol nm^6
C.RbBr.MM  = (215*C_unit)*kj_per_erg*(nm_per_cm^6)*NA; % kJ/mol nm^6

C.RbI.PP  = (59.4*C_unit)*kj_per_erg*(nm_per_cm^6)*NA; % kJ/mol nm^6
C.RbI.PM  = (135*C_unit)*kj_per_erg*(nm_per_cm^6)*NA; % kJ/mol nm^6
C.RbI.MM  = (428*C_unit)*kj_per_erg*(nm_per_cm^6)*NA; % kJ/mol nm^6

C.CsF.PP  = (152*C_unit)*kj_per_erg*(nm_per_cm^6)*NA; % kJ/mol nm^6
C.CsF.PM  = (52*C_unit)*kj_per_erg*(nm_per_cm^6)*NA; % kJ/mol nm^6
C.CsF.MM  = (19.1*C_unit)*kj_per_erg*(nm_per_cm^6)*NA; % kJ/mol nm^6

%% Huggins-Mayer Dipole-Quadrupole Dispersion Parameter D: PM = +-   PP = ++     MM = --
D.LiF.PP  = (0.03*D_unit)*kj_per_erg*(nm_per_cm^8)*NA; % kJ/mol nm^8
D.LiF.PM   = (0.6*D_unit)*kj_per_erg*(nm_per_cm^8)*NA; % kJ/mol nm^8
D.LiF.MM    = (17*D_unit)*kj_per_erg*(nm_per_cm^8)*NA; % kJ/mol nm^8

D.LiCl.PP = (0.03*D_unit)*kj_per_erg*(nm_per_cm^8)*NA; % kJ/mol nm^8
D.LiCl.PM = (2.4*D_unit)*kj_per_erg*(nm_per_cm^8)*NA; % kJ/mol nm^8
D.LiCl.MM   = (223*D_unit)*kj_per_erg*(nm_per_cm^8)*NA; % kJ/mol nm^8

D.LiBr.PP = (0.03*D_unit)*kj_per_erg*(nm_per_cm^8)*NA; % kJ/mol nm^8
D.LiBr.PM = (3.3*D_unit)*kj_per_erg*(nm_per_cm^8)*NA; % kJ/mol nm^8
D.LiBr.MM = (423*D_unit)*kj_per_erg*(nm_per_cm^8)*NA; % kJ/mol nm^8

D.LiI.PP = (0.03*D_unit)*kj_per_erg*(nm_per_cm^8)*NA; % kJ/mol nm^8
D.LiI.PM = (5.3*D_unit)*kj_per_erg*(nm_per_cm^8)*NA; % kJ/mol nm^8
D.LiI.MM = (1060*D_unit)*kj_per_erg*(nm_per_cm^8)*NA; % kJ/mol nm^8

D.NaF.PP  = (0.8*D_unit)*kj_per_erg*(nm_per_cm^8)*NA; % kJ/mol nm^8
D.NaF.PM  = (3.8*D_unit)*kj_per_erg*(nm_per_cm^8)*NA; % kJ/mol nm^8
D.NaF.MM  = (20*D_unit)*kj_per_erg*(nm_per_cm^8)*NA; % kJ/mol nm^8

D.NaCl.PP  = (0.8*D_unit)*kj_per_erg*(nm_per_cm^8)*NA; % kJ/mol nm^8
D.NaCl.PM  = (13.9*D_unit)*kj_per_erg*(nm_per_cm^8)*NA; % kJ/mol nm^8
D.NaCl.MM  = (233*D_unit)*kj_per_erg*(nm_per_cm^8)*NA; % kJ/mol nm^8

D.NaBr.PP  = (0.8*D_unit)*kj_per_erg*(nm_per_cm^8)*NA; % kJ/mol nm^8
D.NaBr.PM  = (19*D_unit)*kj_per_erg*(nm_per_cm^8)*NA; % kJ/mol nm^8
D.NaBr.MM  = (450*D_unit)*kj_per_erg*(nm_per_cm^8)*NA; % kJ/mol nm^8

D.NaI.PP  = (0.8*D_unit)*kj_per_erg*(nm_per_cm^8)*NA; % kJ/mol nm^8
D.NaI.PM  = (31*D_unit)*kj_per_erg*(nm_per_cm^8)*NA; % kJ/mol nm^8
D.NaI.MM  = (1100*D_unit)*kj_per_erg*(nm_per_cm^8)*NA; % kJ/mol nm^8

D.KF.PP  = (24*D_unit)*kj_per_erg*(nm_per_cm^8)*NA; % kJ/mol nm^8
D.KF.PM  = (21*D_unit)*kj_per_erg*(nm_per_cm^8)*NA; % kJ/mol nm^8
D.KF.MM  = (22*D_unit)*kj_per_erg*(nm_per_cm^8)*NA; % kJ/mol nm^8

D.KCl.PP  = (24*D_unit)*kj_per_erg*(nm_per_cm^8)*NA; % kJ/mol nm^8
D.KCl.PM  = (73*D_unit)*kj_per_erg*(nm_per_cm^8)*NA; % kJ/mol nm^8
D.KCl.MM  = (250*D_unit)*kj_per_erg*(nm_per_cm^8)*NA; % kJ/mol nm^8

D.KBr.PP  = (24*D_unit)*kj_per_erg*(nm_per_cm^8)*NA; % kJ/mol nm^8
D.KBr.PM  = (99*D_unit)*kj_per_erg*(nm_per_cm^8)*NA; % kJ/mol nm^8
D.KBr.MM  = (470*D_unit)*kj_per_erg*(nm_per_cm^8)*NA; % kJ/mol nm^8

D.KI.PP  = (24*D_unit)*kj_per_erg*(nm_per_cm^8)*NA; % kJ/mol nm^8
D.KI.PM  = (156*D_unit)*kj_per_erg*(nm_per_cm^8)*NA; % kJ/mol nm^8
D.KI.MM  = (1130*D_unit)*kj_per_erg*(nm_per_cm^8)*NA; % kJ/mol nm^8

D.RbF.PP  = (82*D_unit)*kj_per_erg*(nm_per_cm^8)*NA; % kJ/mol nm^8
D.RbF.PM  = (40*D_unit)*kj_per_erg*(nm_per_cm^8)*NA; % kJ/mol nm^8
D.RbF.MM  = (23*D_unit)*kj_per_erg*(nm_per_cm^8)*NA; % kJ/mol nm^8

D.RbCl.PP  = (82*D_unit)*kj_per_erg*(nm_per_cm^8)*NA; % kJ/mol nm^8
D.RbCl.PM  = (134*D_unit)*kj_per_erg*(nm_per_cm^8)*NA; % kJ/mol nm^8
D.RbCl.MM  = (260*D_unit)*kj_per_erg*(nm_per_cm^8)*NA; % kJ/mol nm^8

D.RbBr.PP  = (82*D_unit)*kj_per_erg*(nm_per_cm^8)*NA; % kJ/mol nm^8
D.RbBr.PM  = (180*D_unit)*kj_per_erg*(nm_per_cm^8)*NA; % kJ/mol nm^8
D.RbBr.MM  = (490*D_unit)*kj_per_erg*(nm_per_cm^8)*NA; % kJ/mol nm^8

D.RbI.PP  = (82*D_unit)*kj_per_erg*(nm_per_cm^8)*NA; % kJ/mol nm^8
D.RbI.PM  = (280*D_unit)*kj_per_erg*(nm_per_cm^8)*NA; % kJ/mol nm^8
D.RbI.MM  = (1200*D_unit)*kj_per_erg*(nm_per_cm^8)*NA; % kJ/mol nm^8

D.CsF.PP  = (278*D_unit)*kj_per_erg*(nm_per_cm^8)*NA; % kJ/mol nm^8
D.CsF.PM  = (78*D_unit)*kj_per_erg*(nm_per_cm^8)*NA; % kJ/mol nm^8
D.CsF.MM  = (23*D_unit)*kj_per_erg*(nm_per_cm^8)*NA; % kJ/mol nm^8

%% TF Repulsive Size Parameter sigma (AKA r+/-): P = +   M = -
% Metals
sigma.Li = 0.816*nm_per_Ang; % nm
sigma.Na = 1.170*nm_per_Ang; % nm
sigma.K  = 1.463*nm_per_Ang; % nm
sigma.Rb = 1.587*nm_per_Ang; % nm
sigma.Cs = 1.720*nm_per_Ang; % nm

% Halides
sigma.F  = 1.179*nm_per_Ang; % nm
sigma.Cl = 1.585*nm_per_Ang; % nm
sigma.Br = 1.716*nm_per_Ang; % nm
sigma.I  = 1.907*nm_per_Ang; % nm

%% TF Parameter: Number of Valence electrons (for Pauling Coefficient Calculation)
% Metals
valence.Li = 2;
valence.Na = 8;
valence.K = 8;
valence.Rb = 8;
valence.Cs = 8;

% Halides
valence.F = 8;
valence.Cl = 8;
valence.Br = 8;
valence.I = 8;

%% TF Hardness Parameter Rho
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

%% Huggins-Mayer potential parameter b (same for all salts)
b = (0.338e-12)*kj_per_erg*NA; % kJ/mol

%% Generate range (r) in nm
r = Startpoint:Spacing:Endpoint;

%% Split Salt Into Component Metal and Halide
[Metal,Halide] = Separate_Metal_Halide(Salt);

%% Calculate Pauling Coefficients beta: PM = +-   PP = ++     MM = --
beta.PP = 1 + 2*q.(Metal)/valence.(Metal); % Unitless
beta.PM = 1 + q.(Metal)/valence.(Metal) + q.(Halide)/valence.(Halide); % Unitless
beta.MM = 1 + 2*q.(Halide)/valence.(Halide); % Unitless

%% Calculate TF Repulsive Exponential Parameter alpha: PM = +-   PP = ++     MM = --
alpha.PP = 1/rho.(Salt); % nm^-1
alpha.PM = 1/rho.(Salt); % nm^-1
alpha.MM = 1/rho.(Salt); % nm^-1

%% Calculate TF Repulsive Scaling Parameter B: PM = +-   PP = ++     MM = --
B.PP = beta.PP*b*exp(2*sigma.(Metal)/rho.(Salt));
B.PM = beta.PM*b*exp((sigma.(Metal) + sigma.(Halide))/rho.(Salt));
B.MM = beta.MM*b*exp(2*sigma.(Halide)/rho.(Salt));

%% Build PES: Plus-Minus

% Plus - Minus total potential
U_MetHal.Total = k_0*(e_c^2)*q.(Metal)*q.(Halide)./(r) ...
    + B.PM*exp(-alpha.PM.*r) ...
    - DS.*C.(Salt).PM./(r.^6) ...
    - DS.*D.(Salt).PM./(r.^8);

% components
U_MetHal.f = 1./r; % Electrostatics function f(r)
U_MetHal.g = - DS.*C.(Salt).PM./(r.^6) - DS.*D.(Salt).PM./(r.^8); % Dispersion g(r)
U_MetHal.h = B.PM*exp(-alpha.PM.*r);% Short range repulsion h(r)

% Plus - Minus total derivative
U_MetHal.dTotal = -k_0*(e_c^2)*q.(Metal)*q.(Halide)./(r.^2) ...
    - alpha.PM*B.PM*exp(-alpha.PM.*r) ...
    + DS.*6.*C.(Salt).PM./(r.^7) ...
    + DS.*8.*D.(Salt).PM./(r.^9);

% components
U_MetHal.df = 1./(r.^2);% Electrostatics function -df(r)/dr
U_MetHal.dg = - DS.*6.*C.(Salt).PM./(r.^7) - DS.*8.*D.(Salt).PM./(r.^9) ; % Dispersion -dg(r)/dr
U_MetHal.dh = alpha.PM*B.PM*exp(-alpha.PM.*r);% Short range repulsion -dh(r)/dr

% remove infinities
U_MetHal = Remove_Infinities(U_MetHal);

%% Build PES: Plus - Plus

% Plus - Plus Total Potential
U_MetMet.Total = k_0*(e_c^2)*q.(Metal).*q.(Metal)./(r) ...
    + B.PP*exp(-alpha.PP.*r) ...
    - DS.*C.(Salt).PP./(r.^6) ...
    - DS.*D.(Salt).PP./(r.^8);

% components
U_MetMet.f = 1./r; % Electrostatics function f(r)
U_MetMet.g = - DS.*C.(Salt).PP./(r.^6) - DS.*D.(Salt).PP./(r.^8); % Dispersion g(r)
U_MetMet.h = B.PP*exp(-alpha.PP.*r); % Short range repulsion

% Plus - Plus Total Derivative
U_MetMet.dTotal = -k_0*(e_c^2)*q.(Metal).*q.(Metal)./(r.^2) ...
    - alpha.PP*B.PP*exp(-alpha.PP.*r)...
    + DS.*6.*C.(Salt).PP./(r.^7) ...
    + DS.*8.*D.(Salt).PP./(r.^9);

% Components
U_MetMet.df = 1./(r.^2); % Electrostatics function f(r)
U_MetMet.dg = - DS.*6.*C.(Salt).PP./(r.^7) - DS.*8.*D.(Salt).PP./(r.^9); % Dispersion g(r)
U_MetMet.dh = alpha.PP*B.PP*exp(-alpha.PP.*r); % Short range repulsion

% Remove infinities
U_MetMet = Remove_Infinities(U_MetMet);

%% Build PES: Minus - Minus
% Minus - Minus Total Potential
U_HalHal.Total = k_0*(e_c^2)*q.(Halide).*q.(Halide)./(r)...
    + B.MM*exp(-alpha.MM.*r) ...
    - DS.*C.(Salt).MM./(r.^6)...
    - DS.*D.(Salt).MM./(r.^8);

% components
U_HalHal.f = 1./r; % Electrostatics function f(r)
U_HalHal.g = - DS.*C.(Salt).MM./(r.^6) - DS.*D.(Salt).MM./(r.^8); % Dispersion g(r)
U_HalHal.h = B.MM*exp(-alpha.MM.*r); % Short range repulsion

% Minus - Minus Total Derivative
U_HalHal.dTotal = -k_0*(e_c^2)*q.(Halide).*q.(Halide)./(r.^2)...
    - alpha.MM*B.MM*exp(-alpha.MM.*r)...
    + DS.*6.*C.(Salt).MM./(r.^7)...
    + DS.*8.*D.(Salt).MM./(r.^9);

% components
U_HalHal.df = 1./(r.^2); % Electrostatics function f(r)
U_HalHal.dg = - DS.*6.*C.(Salt).MM./(r.^7) - DS.*8.*D.(Salt).MM./(r.^9); % Dispersion g(r)
U_HalHal.dh = alpha.MM*B.MM*exp(-alpha.MM.*r); % Short range repulsion

% remove infinities
U_HalHal = Remove_Infinities(U_HalHal);

%% Print
PM = [r ; U_MetHal.f ; U_MetHal.df ; U_MetHal.g ; U_MetHal.dg ; U_MetHal.h ; U_MetHal.dh];
U_PM_out = deblank( sprintf(['%16.10e   %16.10e %16.10e   %16.10e %16.10e   %16.10e %16.10e' newline],PM(:)) );

PP = [r ; U_MetMet.f ; U_MetMet.df ; U_MetMet.g ; U_MetMet.dg ; U_MetMet.h ; U_MetMet.dh];
U_PP_out = deblank( sprintf(['%16.10e   %16.10e %16.10e   %16.10e %16.10e   %16.10e %16.10e' newline],PP(:)) );

MM = [r ; U_HalHal.f ; U_HalHal.df ; U_HalHal.g ; U_HalHal.dg ; U_HalHal.h ; U_HalHal.dh];
U_MM_out = deblank( sprintf(['%16.10e   %16.10e %16.10e   %16.10e %16.10e   %16.10e %16.10e' newline],MM(:)) );

%% PLOT if plotswitch chosen
if plotswitch
    
    % Options
    lw=2;
    fs=25;

    h = cell(1,8);
    hold on

    h{1} = plot(r,U_MetHal.Total,'Color','k','LineWidth',lw,'LineStyle','-');
    h{2} = plot(r,U_MetMet.Total,'Color','r','LineWidth',lw,'LineStyle',':');
    h{3} = plot(r,U_HalHal.Total,'Color','g','LineWidth',lw,'LineStyle','-.');

    title(['Plot of TF Potentials for ' Salt],...
        'Interpreter','latex','fontsize',fs)

    set(gca,'box','on','TickLabelInterpreter','latex');
    set(gca,'XMinorTick','on','YMinorTick','on','FontSize',fs);
    xlabel('r (nm)','fontsize',fs,'Interpreter','latex');
    ylabel('Pair Potential (kJ mol$^{-1}$)','fontsize',fs,'Interpreter','latex');

    ylim([-900 1000]);
    xlim([Startpoint Endpoint]);

    % Blank line
    hline = refline([0 0]);
    hline.Color = 'k';
    hline.LineWidth = lw-1;
    hline.LineStyle = '--';
    leg1 = legend([h{:}],{['TF - ' Salt] ['TF - ' Metal Metal] ['TF - ' Halide Halide]});
    %leg1 = legend([h{:}],{['TF - ' Salt ' V(r)'] ['TF - ' Salt ' V''(r) Analytical'] ['TF - ' Salt ' V''(r) Numerical']});
%     leg1 = legend([h{:}],{['TF - ' Salt ' $V(r)$'] ['TF - ' Salt ' $f(r)$'] ...
%         ['TF - ' Salt ' $g(r)$'] ['TF - ' Salt ' $h(r)$'] ...
%         ['TF - ' Salt ' $\frac{d V (r)}{d r}$'] ['TF - ' Salt ' $\frac{d g (r)}{d r}$'] ...
%         ['TF - ' Salt ' $\frac{d h (r)}{d r}$'] ['TF - ' Salt ' $\frac{d h (r)}{d r}$']});
    leg1.Interpreter = 'latex';
end

end