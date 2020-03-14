% Generates pair potential energy surfaces from 'Startpoint' up to a given
% length input 'Endpoint'. Plotswitch is a logical true or false to
% determine whether to plot the resulting potentials
% Recommended spacing is 0.005 angstroms or 0.0005 nm
% INPUT UNITS MUST BE ALL IN NANOMETERS. OUTPUTS ARE IN NANOMETERS AND kJ/mol

% C6Damping:
% 0 = no (default) damping
% 1 = BJ/rational damping (same as in D3(BJ), damps to a constant. Fairly
% weak damping)
% 2 = Tang Damping (Mid strength damping, damps to zero)
% 3 = MMDRE Damping function (very weak damping)
% 4 = PAMoC Damping function (weak damping)
% 5 = EHFSK Damping function (strong damping)
% 6 = WY damping function (strongest damping)

% C6Damping adds close-range damping

% Parameter sets for C6/C8 coefficients
% 0 = default TF Parameters
% 1 = D3 values
% 2 = Best literature values available
% 3 = D4 with C6 and C8 generated on the fly

% GAdjust are N x 3 arrays of gaussian parameters
% (i , 1) is the Gaussian height of the ith adjustment (may be negative or
% positive)
% (i , 2) is the center point of the ith Gaussian (should be positive)
% (i , 3) is the standard deviation or width (negative and positive values
% are the same)
function [U_MX_out, U_MM_out, U_XX_out] = TFd_Potential_Generator(Startpoint,...
    Endpoint,Spacing,Salt,plotswitch,Scaling_Params,vdw_modifier,RVDW_Cutoff,...
    C6Damping,Paramset,Geometry,Directory,CRDamping,GAdjust_MX,GAdjust_MM,GAdjust_XX)

%% Parameter Scaling 
S_D = Scaling_Params(1); % Dispersion scaling factor
S_R = Scaling_Params(2); % Repulsion scaling factor
S_MMD = S_D*Scaling_Params(5); % Dispersion scaling factor for M-M interaction
S_XXD = S_D*Scaling_Params(6); % Dispersion scaling factor for X-X interaction
S_MXD = S_D*Scaling_Params(7); % Dispersion scaling factor for M-X interaction
S_A = Scaling_Params(8); % Alpha scaling factor

%% Close-Range Sigmoid Damping Parameters
r_d = 0.10; % If close-range dispersion damping is on, this is the value of the sigmoid's midpoint
sb = 150; % sigmoid "steepness" for damping

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
kj_per_erg = 1e-10; % kJ per erg
nm_per_cm = 1e+7; % nm per cm
nm_per_m = 1e+9; % nm per m
NA = 6.0221409e23; % Molecules per mole
e_c = 1.60217662e-19; % Elementary charge in Coulombs
epsilon_0 = (8.854187817620e-12)*1000/(nm_per_m*NA); % Vacuum Permittivity C^2 mol kJ^-1 nm^-1
k_0 = 1/(4*pi*epsilon_0); % Coulomb constant in kJ nm C^-2 mol^-1
C_unit = 1e-60; % erg cm^6
D_unit = 1e-76; % erg cm^8
nm_per_Ang = 0.1; % nm per Angstrom

%% Split Salt Into Component Metal and Halide
[Metal,Halide] = Separate_Metal_Halide(Salt);

if Paramset == 1 % D3 dispersion Parameters
    %% C6 Coefficients
    C.LiF.MM  = 4918.9042107878*(nm_per_Ang^6); % kJ/mol nm^6
    C.LiF.MX  = 990.0857326400*(nm_per_Ang^6); % kJ/mol nm^6
    C.LiF.XX  = 411.2995536808*(nm_per_Ang^6); % kJ/mol nm^6
   
    C.LiCl.MM = 4918.9042107878*(nm_per_Ang^6); % kJ/mol nm^6
    C.LiCl.MX = 4450.9148362278*(nm_per_Ang^6); % kJ/mol nm^6
    C.LiCl.XX = 5211.7103353493*(nm_per_Ang^6); % kJ/mol nm^6
    
    C.LiBr.MM = 4918.9042107878*(nm_per_Ang^6); % kJ/mol nm^6
    C.LiBr.MX = 6238.2901762131*(nm_per_Ang^6); % kJ/mol nm^6
    C.LiBr.XX = 9756.9852141205*(nm_per_Ang^6); % kJ/mol nm^6
    
    C.LiI.MM = 4918.9042107878*(nm_per_Ang^6); % kJ/mol nm^6
    C.LiI.MX = 9306.8796821688*(nm_per_Ang^6); % kJ/mol nm^6
    C.LiI.XX = 20668.4353099615*(nm_per_Ang^6); % kJ/mol nm^6

    C.NaCl.MM  = 10729.4523062025*(nm_per_Ang^6); % kJ/mol nm^6
    C.NaCl.MX  = 6456.1075384258*(nm_per_Ang^6); % kJ/mol nm^6
    C.NaCl.XX  = 5211.710335349*(nm_per_Ang^6); % kJ/mol nm^6
    
    %% C8 coefficients
    D.LiF.MM  = 104130.2216951388*(nm_per_Ang^8); % kJ/mol nm^8
    D.LiF.MX  = 9971.6966373607*(nm_per_Ang^8); % kJ/mol nm^8
    D.LiF.XX  = 1970.7989468780*(nm_per_Ang^8); % kJ/mol nm^8

    D.LiCl.MM = 104130.2216951388*(nm_per_Ang^8); % kJ/mol nm^8
    D.LiCl.MX = 69999.5686578376*(nm_per_Ang^8); % kJ/mol nm^8
    D.LiCl.XX = 60892.5274692634*(nm_per_Ang^8); % kJ/mol nm^8

    D.LiBr.MM = 104130.2216951388*(nm_per_Ang^8); % kJ/mol nm^8
    D.LiBr.MX = 120775.5671985948*(nm_per_Ang^8); % kJ/mol nm^8
    D.LiBr.XX = 172756.4400846956*(nm_per_Ang^8); % kJ/mol nm^8

    D.LiI.MM = 104130.2216951388*(nm_per_Ang^8); % kJ/mol nm^8
    D.LiI.MX = 217169.0376048321*(nm_per_Ang^8); % kJ/mol nm^8
    D.LiI.XX = 531602.2210887696*(nm_per_Ang^8); % kJ/mol nm^8

    D.NaCl.MM  = 390953.8836355778*(nm_per_Ang^8); % kJ/mol nm^8
    D.NaCl.MX  = 133209.9331112395*(nm_per_Ang^8); % kJ/mol nm^8
    D.NaCl.XX  = 60892.5274692634*(nm_per_Ang^8); % kJ/mol nm^8

elseif Paramset == 2 % Best Literature Parameters from doi:10.1080/00268978600102091

    C.LiF.MM  = 4.44*(nm_per_Ang^6); % kJ/mol nm^6
    C.LiF.MX  = 61.46*(nm_per_Ang^6); % kJ/mol nm^6
    C.LiF.XX  = 1103.47*(nm_per_Ang^6); % kJ/mol nm^6

    C.LiCl.MM = 4.44*(nm_per_Ang^6); % kJ/mol nm^6
    C.LiCl.MX = 146.67; % kJ/mol nm^6
    C.LiCl.XX = 8446.11; % kJ/mol nm^6

    C.LiBr.MM = 4.44*(nm_per_Ang^6); % kJ/mol nm^6
    C.LiBr.MX = (2.5*C_unit)*kj_per_erg*(nm_per_cm^6)*NA; % kJ/mol nm^6
    C.LiBr.XX = (185*C_unit)*kj_per_erg*(nm_per_cm^6)*NA; % kJ/mol nm^6

    C.LiI.MM = 4.44*(nm_per_Ang^6); % kJ/mol nm^6
    C.LiI.MX = (3.3*C_unit)*kj_per_erg*(nm_per_cm^6)*NA; % kJ/mol nm^6
    C.LiI.XX = (378*C_unit)*kj_per_erg*(nm_per_cm^6)*NA; % kJ/mol nm^6

    C.NaCl.MM  = 87.29*(nm_per_Ang^6); % kJ/mol nm^6
    C.NaCl.MX  = 716.62*(nm_per_Ang^6); % kJ/mol nm^6
    C.NaCl.XX  = 9195.59*(nm_per_Ang^6); % kJ/mol nm^6

    %% C8 Coefficients
    D.LiF.MM  = 4.05*(nm_per_Ang^8); % kJ/mol nm^8
    D.LiF.MX  = 184.21*(nm_per_Ang^8); % kJ/mol nm^8
    D.LiF.XX  = 5695.74*(nm_per_Ang^8); % kJ/mol nm^8

    D.LiCl.MM = 4.05*(nm_per_Ang^8); % kJ/mol nm^8
    D.LiCl.MX = 760.56*(nm_per_Ang^8); % kJ/mol nm^8
    D.LiCl.XX = 75265.10*(nm_per_Ang^8); % kJ/mol nm^8

    D.LiBr.MM = 4.05*(nm_per_Ang^8); % kJ/mol nm^8
    D.LiBr.MX = (3.3*D_unit)*kj_per_erg*(nm_per_cm^8)*NA; % kJ/mol nm^8
    D.LiBr.XX = (423*D_unit)*kj_per_erg*(nm_per_cm^8)*NA; % kJ/mol nm^8

    D.LiI.MM = 4.05*(nm_per_Ang^8); % kJ/mol nm^8
    D.LiI.MX = (5.3*D_unit)*kj_per_erg*(nm_per_cm^8)*NA; % kJ/mol nm^8
    D.LiI.XX = (1060*D_unit)*kj_per_erg*(nm_per_cm^8)*NA; % kJ/mol nm^8

    D.NaCl.MM  = 198.90*(nm_per_Ang^8); % kJ/mol nm^8
    D.NaCl.MX  = 4388.04*(nm_per_Ang^8); % kJ/mol nm^8
    D.NaCl.XX  = 87244.23*(nm_per_Ang^8); % kJ/mol nm^8
    
elseif Paramset == 3 % D4 on-the-fly parameters
    % Conversion factors for D4
    Bohr_nm = 0.05291772108; % a_0 - > nm
    c6conv = 1e-3/2625.4999/((0.052917726)^6); % J/mol nm^6 - > au (from D3 sourcecode)
    J_kJ = 1e-3; % J - > kJ
    Ha_kJmol = 2625.4999; % Ha - > kJ/mol
    c6units = (1/c6conv)*J_kJ; % au - > kJ/mol nm^6
    c8units = (Ha_kJmol)*(Bohr_nm^8); % au - > kJ/mol nm^8

    % Create vasp format geometry file in given Directory
    N = Geometry.N;   
    geom_txt = ['New Structure' newline '1.0' newline];
    TM = Geometry.Transform*[Geometry.a 0 0; 0 Geometry.b 0; 0 0 Geometry.c];
    for i = 1:3
        for j = 1:3
            geom_txt = [geom_txt pad(num2str(TM(i,j),'%10.10f'),20,'left')];
        end
        geom_txt = [geom_txt newline];
    end
    Nd2 = num2str(N/2);

    geom_txt = [geom_txt pad(Metal,5,'left') pad(Halide,5,'left') newline];
    geom_txt = [geom_txt pad(Nd2,5,'left') pad(Nd2,5,'left') newline 'Direct' newline];
    
    for i = 1:(N/2)
        for j = 1:3
            geom_txt = [geom_txt pad(num2str(Geometry.FC_Metal(i,j),'%10.10f'),16,'left')];
        end
        geom_txt = [geom_txt newline];
    end
    for i = 1:(N/2)
        for j = 1:3
            geom_txt = [geom_txt pad(num2str(Geometry.FC_Halide(i,j),'%10.10f'),16,'left')]; %#ok<*AGROW>
        end
        geom_txt = [geom_txt newline];
    end
    
    % save to file
    filename = [Directory filesep 'tempstruct.vasp'];
    fidPM = fopen(filename,'wt');
    fwrite(fidPM,regexprep(geom_txt,{'\r', '\n\n+'}',{'', '\n'}));
    fclose(fidPM);
    
    % run dftd4 on the structure
    if ispc % for testing
        dftd4 = 'wsl source ~/.bashrc; dftd4 ';
        fnunix = windows2unix(filename);
    elseif isunix
        dftd4 = 'dftd4 ';
        fnunix = filename;
    end
    [ercode,dftd4_out] = system([dftd4 fnunix]);
    
    if ercode ~= 0 %dftd4 failed
        error(['DFTD4 module failed from input: '  dftd4 fnunix])
    end
    
    % Grab C6's from output (in AU)
    C6s = regexp(dftd4_out,'# +Z +covCN +q +C6AA.+\n\n\n','match','ONCE');
    MM_C6 = regexp(C6s,[Metal ' +[-.0-9]+ +[-.0-9]+ +([-.0-9]+) +'],'tokens','once');
    C6_MM = str2double(MM_C6{1});
    XX_C6 = regexp(C6s,[Halide ' +[-.0-9]+ +[-.0-9]+ +([-.0-9]+) +'],'tokens','once');
    C6_XX = str2double(XX_C6{1});
    
    Mol_C6 = regexp(C6s,'Mol\. C6AA.+? + : +([-.0-9]+)','tokens','once');
    C6_Mol = str2double(Mol_C6{1});
    
    % Calculate cross term C6
    C6_MX = (C6_Mol - ((N/2)^2)*C6_MM - ((N/2)^2)*C6_XX)/(2*(N/2)^2);
    
    % Load r2r4 from disc (for calculation of C8)
    loadr2r4 = load('r2r4.mat','r2r4');
    r2r4 = loadr2r4.r2r4;

    % Atomic numbers of atoms in salt
    Z_M = elements('Symbol',Metal,'atomic_number');
    Z_X = elements('Symbol',Halide,'atomic_number');

    %% Calculate C8 for each possible pair of atoms in unit cell, as well as R0AB
    sqrt_Q_M = r2r4(Z_M); % Factor used to calculate C8 for Metal
    sqrt_Q_X = r2r4(Z_X); % Factor used to calculate C8 for Halide

    % C6 coefficients
    C.(Salt).MX = C6_MX*c6units; % in kJ/mol nm^6
    C.(Salt).MM = C6_MM*c6units; % in kJ/mol nm^6
    C.(Salt).XX = C6_XX*c6units; % in kJ/mol nm^6
    
    % C8 coefficients
    D.(Salt).MX = 3.0*(C6_MX)*sqrt_Q_M*sqrt_Q_X*c8units; % in kJ/mol nm^8
    D.(Salt).MM = 3.0*(C6_MM)*sqrt_Q_M*sqrt_Q_M*c8units; % in kJ/mol nm^8
    D.(Salt).XX = 3.0*(C6_XX)*sqrt_Q_X*sqrt_Q_X*c8units; % in kJ/mol nm^8
    
else % Default parameters
    %% Huggins-Mayer Dipole-Dipole Dispersion Parameter C: MX = +-   MM = ++     XX = --
    C.LiF.MM  = (0.073*C_unit)*kj_per_erg*(nm_per_cm^6)*NA; % kJ/mol nm^6
    C.LiF.MX  = (0.8*C_unit)*kj_per_erg*(nm_per_cm^6)*NA; % kJ/mol nm^6
    C.LiF.XX  = (14.5*C_unit)*kj_per_erg*(nm_per_cm^6)*NA; % kJ/mol nm^6

    C.LiCl.MM = (0.073*C_unit)*kj_per_erg*(nm_per_cm^6)*NA; % kJ/mol nm^6
    C.LiCl.MX = (2.0*C_unit)*kj_per_erg*(nm_per_cm^6)*NA; % kJ/mol nm^6
    C.LiCl.XX = (111*C_unit)*kj_per_erg*(nm_per_cm^6)*NA; % kJ/mol nm^6

    C.LiBr.MM = (0.073*C_unit)*kj_per_erg*(nm_per_cm^6)*NA; % kJ/mol nm^6
    C.LiBr.MX = (2.5*C_unit)*kj_per_erg*(nm_per_cm^6)*NA; % kJ/mol nm^6
    C.LiBr.XX = (185*C_unit)*kj_per_erg*(nm_per_cm^6)*NA; % kJ/mol nm^6

    C.LiI.MM = (0.073*C_unit)*kj_per_erg*(nm_per_cm^6)*NA; % kJ/mol nm^6
    C.LiI.MX = (3.3*C_unit)*kj_per_erg*(nm_per_cm^6)*NA; % kJ/mol nm^6
    C.LiI.XX = (378*C_unit)*kj_per_erg*(nm_per_cm^6)*NA; % kJ/mol nm^6

    C.NaF.MM  = (1.68*C_unit)*kj_per_erg*(nm_per_cm^6)*NA; % kJ/mol nm^6
    C.NaF.MX  = (4.5*C_unit)*kj_per_erg*(nm_per_cm^6)*NA; % kJ/mol nm^6
    C.NaF.XX  = (16.5*C_unit)*kj_per_erg*(nm_per_cm^6)*NA; % kJ/mol nm^6

    C.NaCl.MM  = (1.68*C_unit)*kj_per_erg*(nm_per_cm^6)*NA; % kJ/mol nm^6
    C.NaCl.MX  = (11.2*C_unit)*kj_per_erg*(nm_per_cm^6)*NA; % kJ/mol nm^6
    C.NaCl.XX  = (116*C_unit)*kj_per_erg*(nm_per_cm^6)*NA; % kJ/mol nm^6

    C.NaBr.MM  = (1.68*C_unit)*kj_per_erg*(nm_per_cm^6)*NA; % kJ/mol nm^6
    C.NaBr.MX  = (14.0*C_unit)*kj_per_erg*(nm_per_cm^6)*NA; % kJ/mol nm^6
    C.NaBr.XX  = (196*C_unit)*kj_per_erg*(nm_per_cm^6)*NA; % kJ/mol nm^6

    C.NaI.MM  = (1.68*C_unit)*kj_per_erg*(nm_per_cm^6)*NA; % kJ/mol nm^6
    C.NaI.MX  = (19.1*C_unit)*kj_per_erg*(nm_per_cm^6)*NA; % kJ/mol nm^6
    C.NaI.XX  = (392*C_unit)*kj_per_erg*(nm_per_cm^6)*NA; % kJ/mol nm^6

    C.KF.MM  = (24.3*C_unit)*kj_per_erg*(nm_per_cm^6)*NA; % kJ/mol nm^6
    C.KF.MX  = (19.5*C_unit)*kj_per_erg*(nm_per_cm^6)*NA; % kJ/mol nm^6
    C.KF.XX  = (18.6*C_unit)*kj_per_erg*(nm_per_cm^6)*NA; % kJ/mol nm^6

    C.KCl.MM  = (24.3*C_unit)*kj_per_erg*(nm_per_cm^6)*NA; % kJ/mol nm^6
    C.KCl.MX  = (48*C_unit)*kj_per_erg*(nm_per_cm^6)*NA; % kJ/mol nm^6
    C.KCl.XX  = (124.5*C_unit)*kj_per_erg*(nm_per_cm^6)*NA; % kJ/mol nm^6

    C.KBr.MM  = (24.3*C_unit)*kj_per_erg*(nm_per_cm^6)*NA; % kJ/mol nm^6
    C.KBr.MX  = (60*C_unit)*kj_per_erg*(nm_per_cm^6)*NA; % kJ/mol nm^6
    C.KBr.XX  = (206*C_unit)*kj_per_erg*(nm_per_cm^6)*NA; % kJ/mol nm^6

    C.KI.MM  = (24.3*C_unit)*kj_per_erg*(nm_per_cm^6)*NA; % kJ/mol nm^6
    C.KI.MX  = (82*C_unit)*kj_per_erg*(nm_per_cm^6)*NA; % kJ/mol nm^6
    C.KI.XX  = (403*C_unit)*kj_per_erg*(nm_per_cm^6)*NA; % kJ/mol nm^6

    C.RbF.MM  = (59.4*C_unit)*kj_per_erg*(nm_per_cm^6)*NA; % kJ/mol nm^6
    C.RbF.MX  = (31*C_unit)*kj_per_erg*(nm_per_cm^6)*NA; % kJ/mol nm^6
    C.RbF.XX  = (18.9*C_unit)*kj_per_erg*(nm_per_cm^6)*NA; % kJ/mol nm^6

    C.RbCl.MM  = (59.4*C_unit)*kj_per_erg*(nm_per_cm^6)*NA; % kJ/mol nm^6
    C.RbCl.MX  = (79*C_unit)*kj_per_erg*(nm_per_cm^6)*NA; % kJ/mol nm^6
    C.RbCl.XX  = (130*C_unit)*kj_per_erg*(nm_per_cm^6)*NA; % kJ/mol nm^6

    C.RbBr.MM  = (59.4*C_unit)*kj_per_erg*(nm_per_cm^6)*NA; % kJ/mol nm^6
    C.RbBr.MX  = (99*C_unit)*kj_per_erg*(nm_per_cm^6)*NA; % kJ/mol nm^6
    C.RbBr.XX  = (215*C_unit)*kj_per_erg*(nm_per_cm^6)*NA; % kJ/mol nm^6

    C.RbI.MM  = (59.4*C_unit)*kj_per_erg*(nm_per_cm^6)*NA; % kJ/mol nm^6
    C.RbI.MX  = (135*C_unit)*kj_per_erg*(nm_per_cm^6)*NA; % kJ/mol nm^6
    C.RbI.XX  = (428*C_unit)*kj_per_erg*(nm_per_cm^6)*NA; % kJ/mol nm^6

    C.CsF.MM  = (152*C_unit)*kj_per_erg*(nm_per_cm^6)*NA; % kJ/mol nm^6
    C.CsF.MX  = (52*C_unit)*kj_per_erg*(nm_per_cm^6)*NA; % kJ/mol nm^6
    C.CsF.XX  = (19.1*C_unit)*kj_per_erg*(nm_per_cm^6)*NA; % kJ/mol nm^6

    %% Huggins-Mayer Dipole-Quadrupole Dispersion Parameter D: MX = +-   MM = ++     XX = --
    D.LiF.MM  = (0.03*D_unit)*kj_per_erg*(nm_per_cm^8)*NA; % kJ/mol nm^8
    D.LiF.MX   = (0.6*D_unit)*kj_per_erg*(nm_per_cm^8)*NA; % kJ/mol nm^8
    D.LiF.XX    = (17*D_unit)*kj_per_erg*(nm_per_cm^8)*NA; % kJ/mol nm^8

    D.LiCl.MM = (0.03*D_unit)*kj_per_erg*(nm_per_cm^8)*NA; % kJ/mol nm^8
    D.LiCl.MX = (2.4*D_unit)*kj_per_erg*(nm_per_cm^8)*NA; % kJ/mol nm^8
    D.LiCl.XX   = (223*D_unit)*kj_per_erg*(nm_per_cm^8)*NA; % kJ/mol nm^8

    D.LiBr.MM = (0.03*D_unit)*kj_per_erg*(nm_per_cm^8)*NA; % kJ/mol nm^8
    D.LiBr.MX = (3.3*D_unit)*kj_per_erg*(nm_per_cm^8)*NA; % kJ/mol nm^8
    D.LiBr.XX = (423*D_unit)*kj_per_erg*(nm_per_cm^8)*NA; % kJ/mol nm^8

    D.LiI.MM = (0.03*D_unit)*kj_per_erg*(nm_per_cm^8)*NA; % kJ/mol nm^8
    D.LiI.MX = (5.3*D_unit)*kj_per_erg*(nm_per_cm^8)*NA; % kJ/mol nm^8
    D.LiI.XX = (1060*D_unit)*kj_per_erg*(nm_per_cm^8)*NA; % kJ/mol nm^8

    D.NaF.MM  = (0.8*D_unit)*kj_per_erg*(nm_per_cm^8)*NA; % kJ/mol nm^8
    D.NaF.MX  = (3.8*D_unit)*kj_per_erg*(nm_per_cm^8)*NA; % kJ/mol nm^8
    D.NaF.XX  = (20*D_unit)*kj_per_erg*(nm_per_cm^8)*NA; % kJ/mol nm^8

    D.NaCl.MM  = (0.8*D_unit)*kj_per_erg*(nm_per_cm^8)*NA; % kJ/mol nm^8
    D.NaCl.MX  = (13.9*D_unit)*kj_per_erg*(nm_per_cm^8)*NA; % kJ/mol nm^8
    D.NaCl.XX  = (233*D_unit)*kj_per_erg*(nm_per_cm^8)*NA; % kJ/mol nm^8

    D.NaBr.MM  = (0.8*D_unit)*kj_per_erg*(nm_per_cm^8)*NA; % kJ/mol nm^8
    D.NaBr.MX  = (19*D_unit)*kj_per_erg*(nm_per_cm^8)*NA; % kJ/mol nm^8
    D.NaBr.XX  = (450*D_unit)*kj_per_erg*(nm_per_cm^8)*NA; % kJ/mol nm^8

    D.NaI.MM  = (0.8*D_unit)*kj_per_erg*(nm_per_cm^8)*NA; % kJ/mol nm^8
    D.NaI.MX  = (31*D_unit)*kj_per_erg*(nm_per_cm^8)*NA; % kJ/mol nm^8
    D.NaI.XX  = (1100*D_unit)*kj_per_erg*(nm_per_cm^8)*NA; % kJ/mol nm^8

    D.KF.MM  = (24*D_unit)*kj_per_erg*(nm_per_cm^8)*NA; % kJ/mol nm^8
    D.KF.MX  = (21*D_unit)*kj_per_erg*(nm_per_cm^8)*NA; % kJ/mol nm^8
    D.KF.XX  = (22*D_unit)*kj_per_erg*(nm_per_cm^8)*NA; % kJ/mol nm^8

    D.KCl.MM  = (24*D_unit)*kj_per_erg*(nm_per_cm^8)*NA; % kJ/mol nm^8
    D.KCl.MX  = (73*D_unit)*kj_per_erg*(nm_per_cm^8)*NA; % kJ/mol nm^8
    D.KCl.XX  = (250*D_unit)*kj_per_erg*(nm_per_cm^8)*NA; % kJ/mol nm^8

    D.KBr.MM  = (24*D_unit)*kj_per_erg*(nm_per_cm^8)*NA; % kJ/mol nm^8
    D.KBr.MX  = (99*D_unit)*kj_per_erg*(nm_per_cm^8)*NA; % kJ/mol nm^8
    D.KBr.XX  = (470*D_unit)*kj_per_erg*(nm_per_cm^8)*NA; % kJ/mol nm^8

    D.KI.MM  = (24*D_unit)*kj_per_erg*(nm_per_cm^8)*NA; % kJ/mol nm^8
    D.KI.MX  = (156*D_unit)*kj_per_erg*(nm_per_cm^8)*NA; % kJ/mol nm^8
    D.KI.XX  = (1130*D_unit)*kj_per_erg*(nm_per_cm^8)*NA; % kJ/mol nm^8

    D.RbF.MM  = (82*D_unit)*kj_per_erg*(nm_per_cm^8)*NA; % kJ/mol nm^8
    D.RbF.MX  = (40*D_unit)*kj_per_erg*(nm_per_cm^8)*NA; % kJ/mol nm^8
    D.RbF.XX  = (23*D_unit)*kj_per_erg*(nm_per_cm^8)*NA; % kJ/mol nm^8

    D.RbCl.MM  = (82*D_unit)*kj_per_erg*(nm_per_cm^8)*NA; % kJ/mol nm^8
    D.RbCl.MX  = (134*D_unit)*kj_per_erg*(nm_per_cm^8)*NA; % kJ/mol nm^8
    D.RbCl.XX  = (260*D_unit)*kj_per_erg*(nm_per_cm^8)*NA; % kJ/mol nm^8

    D.RbBr.MM  = (82*D_unit)*kj_per_erg*(nm_per_cm^8)*NA; % kJ/mol nm^8
    D.RbBr.MX  = (180*D_unit)*kj_per_erg*(nm_per_cm^8)*NA; % kJ/mol nm^8
    D.RbBr.XX  = (490*D_unit)*kj_per_erg*(nm_per_cm^8)*NA; % kJ/mol nm^8

    D.RbI.MM  = (82*D_unit)*kj_per_erg*(nm_per_cm^8)*NA; % kJ/mol nm^8
    D.RbI.MX  = (280*D_unit)*kj_per_erg*(nm_per_cm^8)*NA; % kJ/mol nm^8
    D.RbI.XX  = (1200*D_unit)*kj_per_erg*(nm_per_cm^8)*NA; % kJ/mol nm^8

    D.CsF.MM  = (278*D_unit)*kj_per_erg*(nm_per_cm^8)*NA; % kJ/mol nm^8
    D.CsF.MX  = (78*D_unit)*kj_per_erg*(nm_per_cm^8)*NA; % kJ/mol nm^8
    D.CsF.XX  = (23*D_unit)*kj_per_erg*(nm_per_cm^8)*NA; % kJ/mol nm^8
end

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

%% Calculate Pauling Coefficients beta: MX = +-   MM = ++     XX = --
beta.MM = 1 + 2*q.(Metal)/valence.(Metal); % Unitless
beta.MX = 1 + q.(Metal)/valence.(Metal) + q.(Halide)/valence.(Halide); % Unitless
beta.XX = 1 + 2*q.(Halide)/valence.(Halide); % Unitless

%% Calculate TF Repulsive Exponential Parameter alpha: MX = +-   MM = ++     XX = --
alpha.MM = S_A/rho.(Salt); % nm^-1
alpha.MX = S_A/rho.(Salt); % nm^-1
alpha.XX = S_A/rho.(Salt); % nm^-1

%% Calculate TF Repulsive Scaling Parameter B: MX = +-   MM = ++     XX = -- (Including scaling)
B.MM = S_R.*beta.MM*b*exp(2*sigma.(Metal)/rho.(Salt));
B.MX = S_R.*beta.MX*b*exp((sigma.(Metal) + sigma.(Halide))/rho.(Salt));
B.XX = S_R.*beta.XX*b*exp(2*sigma.(Halide)/rho.(Salt));

%% Scale Dispersion
C.(Salt).MM = S_MMD.*C.(Salt).MM;
D.(Salt).MM = S_MMD.*D.(Salt).MM;

C.(Salt).XX = S_XXD.*C.(Salt).XX;
D.(Salt).XX = S_XXD.*D.(Salt).XX;

C.(Salt).MX = S_MXD.*C.(Salt).MX;
D.(Salt).MX = S_MXD.*D.(Salt).MX;

%% If Damping at close range
if CRDamping
    f_r = 1./(1 + exp(-sb.*(r - r_d))); % sigmoid damping function
    df_r = (sb.*exp(-sb.*(r - r_d)))./((1 + exp(-sb.*(r - r_d))).^2); % sigmoid damping function derivative
    f_cutoff = 1/(1 + exp(-sb*(RVDW_Cutoff - r_d))); % damping function value at vdw cutoff
else
    f_r = 1; % No damping
    df_r = 0; % Damping derivative is zero
    f_cutoff = 1; % damping function value at vdw cutoff
end

%% No Damping
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
    R0_MX = sqrt(D.(Salt).MX/C.(Salt).MX); % in nm
    R0_MM = sqrt(D.(Salt).MM/C.(Salt).MM); % in nm
    R0_XX = sqrt(D.(Salt).XX/C.(Salt).XX); % in nm
    
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
U_MetHal.Total = f_r.*k_0*(e_c^2)*q.(Metal)*q.(Halide)./(r) ...
    + B.MX*exp(-alpha.MX.*r) ...
    - f6_MX.*C.(Salt).MX./(r.^6) ...
    - f8_MX.*D.(Salt).MX./(r.^8) ...
    + G_r_MX;

% components
U_MetHal.f = 1./r; % Electrostatics function f(r)
U_MetHal.g = - f6_MX.*C.(Salt).MX./(r.^6) - f8_MX.*D.(Salt).MX./(r.^8) + G_r_MX; % Dispersion g(r)
U_MetHal.h = B.MX*exp(-alpha.MX.*r) - k_0*(e_c^2)*q.(Metal)*q.(Halide)./(r) ...
    + f_r.*k_0*(e_c^2)*q.(Metal)*q.(Halide)./(r); % Short range repulsion h(r)

% Plus - Minus total derivative
U_MetHal.dTotal = -f_r.*k_0*(e_c^2)*q.(Metal)*q.(Halide)./(r.^2) ...
    + df_r.*k_0*(e_c^2)*q.(Metal)*q.(Halide)./(r) ...
    - alpha.MX*B.MX*exp(-alpha.MX.*r) ...
    + f6_MX.*6.*C.(Salt).MX./(r.^7) ...
    - df6_MX.*C.(Salt).MX./(r.^6) ...
    + f8_MX.*8.*D.(Salt).MX./(r.^9) ...
    - df8_MX.*D.(Salt).MX./(r.^8) ...
    + dG_r_MX;

% components
U_MetHal.df = 1./(r.^2);% Electrostatics function -df(r)/dr
U_MetHal.dg = - f6_MX.*6.*C.(Salt).MX./(r.^7) + df6_MX.*C.(Salt).MX./(r.^6) ...
    - f8_MX.*8.*D.(Salt).MX./(r.^9) + df8_MX.*D.(Salt).MX./(r.^8) - dG_r_MX; % Dispersion -dg(r)/dr
U_MetHal.dh = alpha.MX*B.MX*exp(-alpha.MX.*r) - k_0*(e_c^2)*q.(Metal)*q.(Halide)./(r.^2) ...
	+ f_r.*k_0*(e_c^2)*q.(Metal)*q.(Halide)./(r.^2) ...
    - df_r.*k_0*(e_c^2)*q.(Metal)*q.(Halide)./(r);% Short range repulsion -dh(r)/dr

if contains(vdw_modifier,'potential-shift','IgnoreCase',true)
    EVDW_Cutoff = B.MX*exp(-alpha.MX*RVDW_Cutoff) ...
        - k_0*(e_c^2)*q.(Metal)*q.(Halide)/(RVDW_Cutoff) ...
        + f_cutoff*k_0*(e_c^2)*q.(Metal)*q.(Halide)/(RVDW_Cutoff) ...
        - f6_cutoff_MX*C.(Salt).MX/(RVDW_Cutoff^6) ...
        - f8_cutoff_MX*D.(Salt).MX/(RVDW_Cutoff^8) ...
        + G_r_MX_Cutoff;
    
    % Shift by the dispersion energy at vdw cutoff radius. only affects one
    % energy component, not derivatives (i.e. forces)
    U_MetHal.Total = U_MetHal.Total - EVDW_Cutoff;
    U_MetHal.g = U_MetHal.g - EVDW_Cutoff;
end

% remove infinities
U_MetHal = Remove_Infinities(U_MetHal);

%% Build PES: Plus - Plus

% Plus - Plus Total Potential
U_MetMet.Total = k_0*(e_c^2)*q.(Metal).*q.(Metal)./(r) ...
    + B.MM*exp(-alpha.MM.*r) ...
    - f6_MM.*C.(Salt).MM./(r.^6) ...
    - f8_MM.*D.(Salt).MM./(r.^8) ...
    + G_r_MM;

% components
U_MetMet.f = 1./r; % Electrostatics function f(r)
U_MetMet.g = - f6_MM.*C.(Salt).MM./(r.^6) - f8_MM.*D.(Salt).MM./(r.^8) + G_r_MM; % Dispersion g(r)
U_MetMet.h = B.MM*exp(-alpha.MM.*r); % Short range repulsion

% Plus - Plus Total Derivative
U_MetMet.dTotal = -k_0*(e_c^2)*q.(Metal).*q.(Metal)./(r.^2) ...
    - alpha.MM*B.MM*exp(-alpha.MM.*r)...
    + f6_MM.*6.*C.(Salt).MM./(r.^7) ...
    - df6_MM.*C.(Salt).MM./(r.^6) ...
    + f8_MM.*8.*D.(Salt).MM./(r.^9) ...
    - df8_MM.*D.(Salt).MM./(r.^8) ...
    + dG_r_MM;

% Components
U_MetMet.df = 1./(r.^2); % Electrostatics function f(r)
U_MetMet.dg = - f6_MM.*6.*C.(Salt).MM./(r.^7) + df6_MM.*C.(Salt).MM./(r.^6)...
    - f8_MM.*8.*D.(Salt).MM./(r.^9) + df8_MM.*D.(Salt).MM./(r.^8) - dG_r_MM; % Dispersion g(r)
U_MetMet.dh = alpha.MM*B.MM*exp(-alpha.MM.*r); % Short range repulsion

if contains(vdw_modifier,'potential-shift','IgnoreCase',true)
    EVDW_Cutoff = B.MM*exp(-alpha.MM.*RVDW_Cutoff) ...
    - f6_cutoff_MM.*C.(Salt).MM./(RVDW_Cutoff.^6) ...
    - f8_cutoff_MM.*D.(Salt).MM./(RVDW_Cutoff.^8) ...
    + G_r_MM_Cutoff;
    
    % Shift by the dispersion energy at vdw cutoff radius. only affects one
    % energy component, not derivatives (i.e. forces)
    U_MetMet.Total = U_MetMet.Total - EVDW_Cutoff;
    U_MetMet.g = U_MetMet.g - EVDW_Cutoff;
end

% Remove infinities
U_MetMet = Remove_Infinities(U_MetMet);

%% Build PES: Minus - Minus
% Minus - Minus Total Potential
U_HalHal.Total = k_0*(e_c^2)*q.(Halide).*q.(Halide)./(r)...
    + B.XX*exp(-alpha.XX.*r) ...
    - f6_XX.*C.(Salt).XX./(r.^6)...
    - f8_XX.*D.(Salt).XX./(r.^8) ...
    + G_r_XX;

% components
U_HalHal.f = 1./r; % Electrostatics function f(r)
U_HalHal.g = - f6_XX.*C.(Salt).XX./(r.^6) - f8_XX.*D.(Salt).XX./(r.^8) + G_r_XX; % Dispersion g(r)
U_HalHal.h = B.XX*exp(-alpha.XX.*r); % Short range repulsion

% Minus - Minus Total Derivative
U_HalHal.dTotal = -k_0*(e_c^2)*q.(Halide).*q.(Halide)./(r.^2)...
    - alpha.XX*B.XX*exp(-alpha.XX.*r)...
    + f6_XX.*6.*C.(Salt).XX./(r.^7) ...
    - df6_XX.*C.(Salt).XX./(r.^6) ...
    + f8_XX.*8.*D.(Salt).XX./(r.^9) ...
    - df8_XX.*D.(Salt).XX./(r.^8) ...
    + dG_r_XX;

% components
U_HalHal.df = 1./(r.^2); % Electrostatics function f(r)
U_HalHal.dg = - f6_XX.*6.*C.(Salt).XX./(r.^7) + df6_XX.*C.(Salt).XX./(r.^6) ...
    - f8_XX.*8.*D.(Salt).XX./(r.^9) + df8_XX.*D.(Salt).XX./(r.^8) - dG_r_XX; % Dispersion g(r)
U_HalHal.dh = alpha.XX*B.XX*exp(-alpha.XX.*r); % Short range repulsion

if contains(vdw_modifier,'potential-shift','IgnoreCase',true)
    EVDW_Cutoff = B.XX*exp(-alpha.XX.*RVDW_Cutoff) ...
    - f6_cutoff_XX.*C.(Salt).XX./(RVDW_Cutoff.^6)...
    - f8_cutoff_XX.*D.(Salt).XX./(RVDW_Cutoff.^8) ...
    + G_r_XX_Cutoff;
    
    % Shift by the dispersion energy at vdw cutoff radius. only affects one
    % energy component, not derivatives (i.e. forces)
    U_HalHal.Total = U_HalHal.Total - EVDW_Cutoff;
    U_HalHal.g = U_HalHal.g - EVDW_Cutoff;
end

% remove infinities
U_HalHal = Remove_Infinities(U_HalHal);

%% Print to table
MX = [r ; U_MetHal.f ; U_MetHal.df ; U_MetHal.g ; U_MetHal.dg ; U_MetHal.h ; U_MetHal.dh];
U_MX_out = deblank( sprintf(['%16.10e   %16.10e %16.10e   %16.10e %16.10e   %16.10e %16.10e' newline],MX(:)) );

MM = [r ; U_MetMet.f ; U_MetMet.df ; U_MetMet.g ; U_MetMet.dg ; U_MetMet.h ; U_MetMet.dh];
U_MM_out = deblank( sprintf(['%16.10e   %16.10e %16.10e   %16.10e %16.10e   %16.10e %16.10e' newline],MM(:)) );

XX = [r ; U_HalHal.f ; U_HalHal.df ; U_HalHal.g ; U_HalHal.dg ; U_HalHal.h ; U_HalHal.dh];
U_XX_out = deblank( sprintf(['%16.10e   %16.10e %16.10e   %16.10e %16.10e   %16.10e %16.10e' newline],XX(:)) );

%% PLOT if plotswitch chosen
if plotswitch
    figure;
    
    % Options
    lw=2;
    fs=25;

    h = cell(1,8);
    hold on

    h{1} = plot(r,U_MetHal.g+U_MetHal.h,'Color','b','LineWidth',lw,'LineStyle','-');
    h{2} = plot(r,U_MetMet.g+U_MetMet.h,'Color','r','LineWidth',lw,'LineStyle','-');
    h{3} = plot(r,U_HalHal.g+U_HalHal.h,'Color','g','LineWidth',lw,'LineStyle','-');
    
%     dr = r(2) - r(1);
%     r_d = r(2:end)- dr/2;
%     
% %     h{1} = plot(r,-U_MetHal.dg,'Color','b','LineWidth',lw,'LineStyle',':');
% %     h{2} = plot(r,-U_MetMet.dg,'Color','r','LineWidth',lw,'LineStyle',':');
% %     h{3} = plot(r,-U_HalHal.dg,'Color','g','LineWidth',lw,'LineStyle',':');
%     
%     h{4} = plot(r_d,diff(U_MetHal.g)/dr,'Color','k','LineWidth',lw,'LineStyle',':');
%     h{5} = plot(r_d,diff(U_MetMet.g)/dr,'Color','k','LineWidth',lw,'LineStyle',':');
%     h{6} = plot(r_d,diff(U_HalHal.g)/dr,'Color','k','LineWidth',lw,'LineStyle',':');
    

    title(['Plot of TF Potentials for ' Salt],...
        'Interpreter','latex','fontsize',fs)

    set(gca,'box','on','TickLabelInterpreter','latex');
    set(gca,'XMinorTick','on','YMinorTick','on','FontSize',fs);
    xlabel('r (nm)','fontsize',fs,'Interpreter','latex');
    ylabel('Pair Potential (kJ mol$^{-1}$)','fontsize',fs,'Interpreter','latex');

    ylim([-100 500]);
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