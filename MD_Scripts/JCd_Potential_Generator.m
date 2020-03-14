% Jung-Cheatham model parameters adapted for three water models with
% Lorenz-Berthelot combining rules
% Use Watermodel = 'SPC/E', 'TIP3P', or 'TIP4PEW'
% Generates pair potential energy surfaces from 'Startpoint' up to a given
% length input 'Endpoint'. Plotswitch is a logical true or false to
% determine whether to plot the resulting potentials
% Recommended spacing is 0.005 angstroms or 0.0005 nm
% INPUT UNITS MUST BE ALL IN NANOMETERS. OUTPUTS ARE IN NANOMETER AND
% kJ/mol

% C6Damping:
% 0 = no (default) damping. This is default of JC model.
% 1 = BJ/rational damping (same as in D3(BJ))
% 2 = Tang Damping (Essentially removes dispersion in JC)
% 3 = MMDRE Damping function (Very weak damping, damps mainly at mid range)
% 4 = PAMoC Damping function (fairly weak damping, damps mainly at mid range)
% 5 = EHFSK Damping function (strong damping)
% 6 = WY damping function (strongest damping)

% GAdjust are N x 3 arrays of gaussian parameters
% (i , 1) is the Gaussian height of the ith adjustment (may be negative or
% positive)
% (i , 2) is the center point of the ith Gaussian (should be positive)
% (i , 3) is the standard deviation or width (negative and positive values
% are the same)

% CRDamping = Close Range Damping
function [U_MX_out, U_MM_out, U_XX_out] = JCd_Potential_Generator(Startpoint,...
    Endpoint,Spacing,Salt,Watermodel,plotswitch,Scaling_Params,vdw_modifier,RVDW_Cutoff,C6Damping,CRDamping,...
    GAdjust_MX,GAdjust_MM,GAdjust_XX)


S_D = Scaling_Params(1);
S_R = Scaling_Params(2);
S_E = Scaling_Params(3);
S_S = Scaling_Params(4);
S_MMD = Scaling_Params(5);
S_XXD = Scaling_Params(6);
S_MXD = Scaling_Params(7);
% Ignore 8th scaling parameter (only applies to TF model)

%% Close-Range Sigmoid Damping Parameters
r_d = 0.10; % If close-range dispersion damping is on, this is the value of the sigmoid's midpoint in nm
b = 150; % sigmoid "steepness" for damping

%% Gaussian adjustments
if exist('GAdjust_MM','var')
    G_a_MM = GAdjust_MM(:,1);
    G_b_MM = GAdjust_MM(:,2);
    G_c_MM = GAdjust_MM(:,3);
else
    G_a_MM = 0;
    G_b_MM = 0;
    G_c_MM = 1;
end

if exist('GAdjust_XX','var')
    G_a_XX = GAdjust_XX(:,1);
    G_b_XX = GAdjust_XX(:,2);
    G_c_XX = GAdjust_XX(:,3);
else
    G_a_XX = 0;
    G_b_XX = 0;
    G_c_XX = 1;
end

if exist('GAdjust_MX','var')
    G_a_MX = GAdjust_MX(:,1);
    G_b_MX = GAdjust_MX(:,2);
    G_c_MX = GAdjust_MX(:,3);
else
    G_a_MX = 0;
    G_b_MX = 0;
    G_c_MX = 1;
end

%% Conversion factors and fundamental constants
nm_per_m = 1e+9; % nm per m
NA = 6.0221409e23; % Molecules per mole
e_c = 1.60217662e-19; % Elementary charge in Coulombs
epsilon_0 = (8.854187817620e-12)*1000/(nm_per_m*NA); % Vacuum Permittivity C^2 mol kJ^-1 nm^-1
k_0 = 1/(4*pi*epsilon_0); % Coulomb constant in kJ nm C^-2 mol^-1
nm_per_Ang = 0.1; % nm per Angstrom
kJ_per_kcal = 4.184; % kj per kcal

[Metal,Halide] = Separate_Metal_Halide(Salt);

%% JC Ion Parameters in SPC/E water
if strcmp(Watermodel,'SPC/E')
    Param.Li.sigma = (0.791*nm_per_Ang)*(2^(5/6)); % nm
    Param.Li.epsilon = (0.3367344)*kJ_per_kcal; % kJ/mol

    Param.Na.sigma = (1.212*nm_per_Ang)*(2^(5/6)); % nm
    Param.Na.epsilon = (0.3526418)*kJ_per_kcal; % kJ/mol

    Param.K.sigma = (1.593*nm_per_Ang)*(2^(5/6)); % nm
    Param.K.epsilon = (0.4297054)*kJ_per_kcal; % kJ/mol

    Param.Rb.sigma = (1.737*nm_per_Ang)*(2^(5/6)); % nm
    Param.Rb.epsilon = (0.4451036)*kJ_per_kcal; % kJ/mol

    Param.Cs.sigma = (2.021*nm_per_Ang)*(2^(5/6)); % nm
    Param.Cs.epsilon = (0.0898565)*kJ_per_kcal; % kJ/mol

    Param.F.sigma = (2.257*nm_per_Ang)*(2^(5/6)); % nm
    Param.F.epsilon = (0.0074005)*kJ_per_kcal; % kJ/mol

    Param.Cl.sigma = (2.711*nm_per_Ang)*(2^(5/6)); % nm
    Param.Cl.epsilon = (0.0127850)*kJ_per_kcal; % kJ/mol

    Param.Br.sigma = (2.751*nm_per_Ang)*(2^(5/6)); % nm
    Param.Br.epsilon = (0.0269586)*kJ_per_kcal; % kJ/mol

    Param.I.sigma = (2.919*nm_per_Ang)*(2^(5/6)); % nm
    Param.I.epsilon = (0.0427845)*kJ_per_kcal; % kJ/mol

elseif strcmp(Watermodel,'TIP3P')
    %% JC Ion Parameters in TIP3P water
    Param.Li.sigma = (1.025*nm_per_Ang)*(2^(5/6)); % nm
    Param.Li.epsilon = (0.0279896)*kJ_per_kcal; % kJ/mol

    Param.Na.sigma = (1.369*nm_per_Ang)*(2^(5/6)); % nm
    Param.Na.epsilon = (0.0874393)*kJ_per_kcal; % kJ/mol

    Param.K.sigma = (1.705*nm_per_Ang)*(2^(5/6)); % nm
    Param.K.epsilon = (0.1936829)*kJ_per_kcal; % kJ/mol

    Param.Rb.sigma = (1.813*nm_per_Ang)*(2^(5/6)); % nm
    Param.Rb.epsilon = (0.3278219)*kJ_per_kcal; % kJ/mol

    Param.Cs.sigma = (1.976*nm_per_Ang)*(2^(5/6)); % nm
    Param.Cs.epsilon = (0.4065394)*kJ_per_kcal; % kJ/mol

    Param.F.sigma = (2.303*nm_per_Ang)*(2^(5/6)); % nm
    Param.F.epsilon = (0.0033640)*kJ_per_kcal; % kJ/mol

    Param.Cl.sigma = (2.513*nm_per_Ang)*(2^(5/6)); % nm
    Param.Cl.epsilon = (0.0355910)*kJ_per_kcal; % kJ/mol

    Param.Br.sigma = (2.608*nm_per_Ang)*(2^(5/6)); % nm
    Param.Br.epsilon = (0.0586554)*kJ_per_kcal; % kJ/mol

    Param.I.sigma = (2.860*nm_per_Ang)*(2^(5/6)); % nm
    Param.I.epsilon = (0.0536816)*kJ_per_kcal; % kJ/mol

elseif strcmp(Watermodel,'TIP4PEW')
    %% JC Ion Parameters in TIP4P water
    Param.Li.sigma = (0.808*nm_per_Ang)*(2^(5/6)); % nm
    Param.Li.epsilon = (0.1039884)*kJ_per_kcal; % kJ/mol

    Param.Na.sigma = (1.226*nm_per_Ang)*(2^(5/6)); % nm
    Param.Na.epsilon = (0.1684375)*kJ_per_kcal; % kJ/mol

    Param.K.sigma = (1.590*nm_per_Ang)*(2^(5/6)); % nm
    Param.K.epsilon = (0.2794651)*kJ_per_kcal; % kJ/mol

    Param.Rb.sigma = (1.709*nm_per_Ang)*(2^(5/6)); % nm
    Param.Rb.epsilon = (0.4331494)*kJ_per_kcal; % kJ/mol

    Param.Cs.sigma = (1.888*nm_per_Ang)*(2^(5/6)); % nm
    Param.Cs.epsilon = (0.3944318)*kJ_per_kcal; % kJ/mol

    Param.F.sigma = (2.538*nm_per_Ang)*(2^(5/6)); % nm
    Param.F.epsilon = (0.0015752)*kJ_per_kcal; % kJ/mol

    Param.Cl.sigma = (2.760*nm_per_Ang)*(2^(5/6)); % nm
    Param.Cl.epsilon = (0.0116615)*kJ_per_kcal; % kJ/mol

    Param.Br.sigma = (2.768*nm_per_Ang)*(2^(5/6)); % nm
    Param.Br.epsilon = (0.0303773)*kJ_per_kcal; % kJ/mol

    Param.I.sigma = (2.952*nm_per_Ang)*(2^(5/6)); % nm
    Param.I.epsilon = (0.0417082)*kJ_per_kcal; % kJ/mol
else
    error(['Unknown Water Model: "' Watermodel ...
        '". Please choose one of "SPC/E", "TIP3P", or "TIP4PEW".'])
end

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


%% Calculate parameters of interest for LJ potential
sigma_MM = S_S*Param.(Metal).sigma;
sigma_XX = S_S*Param.(Halide).sigma;
sigma_MX = ( sigma_MM + sigma_XX )/2;

epsilon_MM = S_E*Param.(Metal).epsilon;
epsilon_XX = S_E*Param.(Halide).epsilon;
epsilon_MX = sqrt(epsilon_MM*epsilon_XX);

% Change parameteters into A/r12 - B/r6 format
A_MM = S_R*4*epsilon_MM*(sigma_MM^12);
B_MM = S_MMD*S_D*4*epsilon_MM*(sigma_MM^6);

A_XX = S_R*4*epsilon_XX*(sigma_XX^12);
B_XX = S_XXD*S_D*4*epsilon_XX*(sigma_XX^6);

A_MX = S_R*4*epsilon_MX*(sigma_MX^12);
B_MX = S_MXD*S_D*4*epsilon_MX*(sigma_MX^6);

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

%% No Damping
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

    % use the TF hardness parameters
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

%% Modify potential with Gaussian Adjustments
G_r_MM = zeros(1,length(r));
dG_r_MM = zeros(1,length(r));
G_r_MM_Cutoff = 0;
for i = 1:length(G_a_MM)
    G_r_MM = G_r_MM + G_a_MM(i).*exp((-(r - G_b_MM(i)).^2)./(2.*(G_c_MM(i).^2)));
    G_r_MM_Cutoff = G_r_MM_Cutoff + G_a_MM(i)*exp((-(RVDW_Cutoff - G_b_MM(i))^2)/(2*(G_c_MM(i)^2)));
    dG_r_MM = dG_r_MM - (G_a_MM(i).*(r - G_b_MM(i))).*(exp((-(r - G_b_MM(i)).^2)./(2.*(G_c_MM(i).^2))))/(G_c_MM(i).^2);
end

G_r_XX = zeros(1,length(r));
dG_r_XX = zeros(1,length(r));
G_r_XX_Cutoff = 0;
for i = 1:length(G_a_XX)
    G_r_XX = G_r_XX + G_a_XX(i).*exp((-(r - G_b_XX(i)).^2)./(2.*(G_c_XX(i).^2)));
	G_r_XX_Cutoff = G_r_XX_Cutoff + G_a_XX(i)*exp((-(RVDW_Cutoff - G_b_XX(i))^2)/(2*(G_c_XX(i)^2)));
    dG_r_XX = dG_r_XX - (G_a_XX(i).*(r - G_b_XX(i))).*(exp((-(r - G_b_XX(i)).^2)./(2.*(G_c_XX(i).^2))))/(G_c_XX(i).^2);
end

G_r_MX = zeros(1,length(r));
dG_r_MX = zeros(1,length(r));
G_r_MX_Cutoff = 0;
for i = 1:length(G_a_MX)
    G_r_MX = G_r_MX + G_a_MX(i).*exp((-(r - G_b_MX(i)).^2)./(2.*(G_c_MX(i).^2)));
    G_r_MX_Cutoff = G_r_MX_Cutoff + G_a_MX(i)*exp((-(RVDW_Cutoff - G_b_MX(i))^2)/(2*(G_c_MX(i)^2)));
    dG_r_MX = dG_r_MX - (G_a_MX(i).*(r - G_b_MX(i))).*(exp((-(r - G_b_MX(i)).^2)./(2.*(G_c_MX(i).^2))))/(G_c_MX(i).^2);
end

%% Build PES: Plus-Minus

% Plus - Minus total potential
U_MX.Total = f_r.*k_0*(e_c^2)*q.(Metal)*q.(Halide)./(r) ...
    + A_MX./(r.^12) ...
    - f6_MX.*B_MX./(r.^6) ...
    + G_r_MX;

% components
U_MX.f = 1./r; % Electrostatics function f(r)
U_MX.g = -f6_MX.*B_MX./(r.^6) + G_r_MX; % Dispersion g(r)
U_MX.h = A_MX./(r.^12) - k_0*(e_c^2)*q.(Metal)*q.(Halide)./(r) + ...
    f_r.*k_0*(e_c^2)*q.(Metal)*q.(Halide)./(r); % Short range repulsion

% Plus - Minus total derivative
U_MX.dTotal = -f_r.*k_0*(e_c^2)*q.(Metal)*q.(Halide)./(r.^2) ...
    + df_r.*k_0*(e_c^2)*q.(Metal)*q.(Halide)./(r)...
    - A_MX.*12./(r.^13) ...
    + f6_MX.*B_MX.*6./(r.^7) ...
    - df6_MX.*B_MX./(r.^6) ...
    + dG_r_MX;

% components
U_MX.df = 1./(r.^2);% Electrostatics function
U_MX.dg = - f6_MX.*B_MX.*6./(r.^7) + df6_MX.*B_MX./(r.^6) - dG_r_MX; % Dispersion -dg(r)/dr
U_MX.dh = + A_MX.*12./(r.^13) - k_0*(e_c^2)*q.(Metal)*q.(Halide)./(r.^2) ...
    + f_r.*k_0*(e_c^2)*q.(Metal)*q.(Halide)./(r.^2) ...
    - df_r.*k_0*(e_c^2)*q.(Metal)*q.(Halide)./(r);% Short range repulsion

if contains(vdw_modifier,'potential-shift','IgnoreCase',true)
    EVDW_Cutoff = A_MX./(RVDW_Cutoff.^12) ...
        - k_0*(e_c^2)*q.(Metal)*q.(Halide)/(RVDW_Cutoff) ...
        + f_cutoff*k_0*(e_c^2)*q.(Metal)*q.(Halide)/(RVDW_Cutoff) ...
        - f6_cutoff_MX.*B_MX./(RVDW_Cutoff.^6) ...
        + G_r_MX_Cutoff;
    
    % Shift by the dispersion energy at vdw cutoff radius. only affects one
    % energy component, not derivatives (i.e. forces)
    U_MX.Total = U_MX.Total - EVDW_Cutoff;
    U_MX.g = U_MX.g - EVDW_Cutoff;
end

% remove infinities
U_MX = Remove_Infinities(U_MX);

%% Build PES: Plus - Plus

% Plus - Plus total potential
U_MM.Total = k_0*(e_c^2)*q.(Metal)*q.(Metal)./(r) ...
    + A_MM./(r.^12) ...
    - f6_MM.*B_MM./(r.^6) ...
    + G_r_MM;

% components
U_MM.f = 1./r; % Electrostatics function f(r)
U_MM.g = -f6_MM.*B_MM./(r.^6) + G_r_MM; % Dispersion g(r)
U_MM.h = A_MM./(r.^12);% Short range repulsion

% Plus - Plus total derivative
U_MM.dTotal = -k_0*(e_c^2)*q.(Metal)*q.(Metal)./(r.^2) ...
    - A_MM.*12./(r.^13) ...
    + f6_MM.*B_MM.*6./(r.^7) ...
    - df6_MM.*B_MM./(r.^6) ...
    + dG_r_MM;

% components
U_MM.df = 1./(r.^2);% Electrostatics function
U_MM.dg = - f6_MM.*B_MM.*6./(r.^7) + df6_MM.*B_MM./(r.^6) - dG_r_MM; % Dispersion
U_MM.dh = + A_MM.*12./(r.^13);% Short range repulsion

if contains(vdw_modifier,'potential-shift','IgnoreCase',true)
    EVDW_Cutoff = A_MM./(RVDW_Cutoff.^12) ...
    - f6_cutoff_MM.*B_MM./(RVDW_Cutoff.^6) ...
    + G_r_MM_Cutoff;
    
    % Shift by the dispersion energy at vdw cutoff radius. only affects one
    % energy component, not derivatives (i.e. forces)
    U_MM.Total = U_MM.Total - EVDW_Cutoff;
    U_MM.g = U_MM.g - EVDW_Cutoff;
end

% remove infinities
U_MM = Remove_Infinities(U_MM);

%% Build PES: Minus - Minus

% Minus - Minus total potential
U_XX.Total = k_0*(e_c^2)*q.(Halide)*q.(Halide)./(r) ...
    + A_XX./(r.^12) ...
    - f6_XX.*B_XX./(r.^6) ...
    + G_r_XX;

% components
U_XX.f = 1./r; % Electrostatics function f(r)
U_XX.g = -f6_XX.*B_XX./(r.^6) + G_r_XX; % Dispersion g(r)
U_XX.h = A_XX./(r.^12);% Short range repulsion

% Minus - Minus total derivative
U_XX.dTotal = -k_0*(e_c^2)*q.(Halide)*q.(Halide)./(r.^2) ...
    - A_XX.*12./(r.^13) ...
    + f6_XX.*B_XX.*6./(r.^7) ...
    - df6_XX.*B_XX./(r.^6) ...
    + dG_r_XX;

% components
U_XX.df = 1./(r.^2);% Electrostatics function
U_XX.dg = - f6_XX.*B_XX.*6./(r.^7) + df6_XX.*B_XX./(r.^6) - dG_r_XX; % Dispersion
U_XX.dh = + A_XX.*12./(r.^13);% Short range repulsion

if contains(vdw_modifier,'potential-shift','IgnoreCase',true)
    EVDW_Cutoff = A_XX./(RVDW_Cutoff.^12) ...
    - f6_cutoff_XX.*B_XX./(RVDW_Cutoff.^6) ...
    + G_r_XX_Cutoff;
    
    % Shift by the dispersion energy at vdw cutoff radius. only affects one
    % energy component, not derivatives (i.e. forces)
    U_XX.Total = U_XX.Total - EVDW_Cutoff;
    U_XX.g = U_XX.g - EVDW_Cutoff;
end

% remove infinities
U_XX = Remove_Infinities(U_XX);

%% Print
MX = [r ; U_MX.f ; U_MX.df ; U_MX.g ; U_MX.dg ; U_MX.h ; U_MX.dh];
U_MX_out = deblank( sprintf(['%16.10e   %16.10e %16.10e   %16.10e %16.10e   %16.10e %16.10e' newline],MX(:)) );

MM = [r ; U_MM.f ; U_MM.df ; U_MM.g ; U_MM.dg ; U_MM.h ; U_MM.dh];
U_MM_out = deblank( sprintf(['%16.10e   %16.10e %16.10e   %16.10e %16.10e   %16.10e %16.10e' newline],MM(:)) );

XX = [r ; U_XX.f ; U_XX.df ; U_XX.g ; U_XX.dg ; U_XX.h ; U_XX.dh];
U_XX_out = deblank( sprintf(['%16.10e   %16.10e %16.10e   %16.10e %16.10e   %16.10e %16.10e' newline],XX(:)) );


%% PLOT if plotswitch chosen
if plotswitch
    figure;
    % Options
    lw=2;
    fs=25;

    h = cell(1,8);
    hold on

    h{1} = plot(r,U_MX.Total,'Color','r','LineWidth',lw,'LineStyle','-');
    h{2} = plot(r,U_MM.Total,'Color','b','LineWidth',lw,'Linestyle','-');
    h{3} = plot(r,U_XX.Total,'Color','g','LineWidth',lw,'Linestyle','-');

%     GOI = gradient(U_MX.Total,Spacing);

    title(['Plot of JC Potentials for ' Salt],...
        'Interpreter','latex','fontsize',fs)

    set(gca,'box','on','TickLabelInterpreter','latex');
    set(gca,'XMinorTick','on','YMinorTick','on','FontSize',fs);
    xlabel('r (nm)','fontsize',fs,'Interpreter','latex');
    ylabel('Pair Potential (kJ mol$^{-1}$)','fontsize',fs,'Interpreter','latex');

     ylim([-1000 1000]);
%    ylim([-100 50]);
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