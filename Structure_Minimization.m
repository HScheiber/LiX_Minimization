%% INFO ABOUT INPUTS
% Models is a string variable. One of 'JC' or 'TF'. Sets the mathematical form of the model to use.
% Salts is either a string or cell array. Can be any of: 'LiF' 'LiCl' 'LiBr' 'LiI' 'NaCl'
% Structures is either a string or cell array. Can be any of: 'Rocksalt' 'Wurtzite' 'Sphalerite' 'CsCl' 'NiAs' 'BetaBeO' 'FiveFive'
% OptPos is a Boolean. If true, then the optimization algorithm optimizes
% for lattice parameters AND fractional coordinates. If false, the
% algorithm only optimizes for lattice parameters.
%% Parameters is and N x M matrix of floats, N and M depend on the chosen model. 
% For JC model: Parameters are contained in either a [2 x 2] OR [2 x 3] array.
% For TF model: Parameters are contained in a [4 x 3] array.
%
% Below are descriptions of the Parameter matrix options. M refers to the
% Metal ion (i.e. Li or Na) and X refers to the Halide ion (i.e. F, Cl, Br,
% or I)
%
% JC 2 x 2 matrix:   sigma_M      sigma_X                (units: nm)
%                    epsilon_M    epsilon_X              (units: kJ/mol)
%
% When using a 2 x 2 array with the JC model, 
% the Lorenz-Berthelot mixing rules are assumed.
% 
%
% JC 2 x 3 array:   sigma_MM    sigma_XX    sigma_MX         (units: nm)
%                   epsilon_MM  epsilon_XX  epsilon_MX       (units: kJ/mol)
%
% When using a 2 x 3 array with the JC model, mixing rules are defined by the input.
%
%
% TF 4 x 3 array:   alpha_MM    alpha_XX    alpha_MX   (units: nm^-1)
%                   B_MM        B_XX        B_MX       (units: kJ/mol)
%                   C_MM        C_XX        C_MX       (units: (kJ nm^6)/mol)
%                   D_MM        D_XX        D_MX       (units: (kJ nm^8)/mol)
%
% No combining rules are defined for the TF model. Note that in original TF
% model, all alpha parameters (which correspond to repulsive wall steepness) are set equal for each salt.
function Output_Array = Structure_Minimization(Salt,Structure,Model,Parameters,OptPos)
%% Structure Settings
Parallel_Mode = true; % If set to true, this will pin the GROMACS process to a single thread
Data_Types = 1; % Allowed data types for automatic search of initial conditions (0 = not optimized, 1 = cell optimized, 2 = full optimized, 3 = atom optimized only)
Continue_if_no_IC = true; % When true, uses input initial conditions if none are found in data. When false, does not attempt calculation if no IC are in data.
Find_Min_Params = true; % When true, finds lowest energy parameters for IC based on Data_Types. When false, uses input IC
Find_Similar_Params = false; % When true, finds lowest energy parameters for IC if possible, but if no data is available, also looks for the same IC with non-scaled model
N_atoms = 10000; % Minimum number of atoms to include in super cell
emtol = 1e-3; %1e-3 [kJ mol-1 nm-1] The position minimization is converged (per cycle) when the maximum force is smaller than this value
MaxCycles = 1000; % Maximum number of optimization cycles.
Maintain_Symmetry = true; % If true, maintains the space group symmetry of the unit cell
Energy_Tol = 1e-4; % kJ/mol. Convergence reached when change in energy between cycles is less than this value
Gradient_Tol_RMS = 0.58959; %kJ/(mol Ang) Convergence criteria 2: must have RMS cell parameter gradient less than this value for convergence
Gradient_Tol_Max = 0.88439; %kJ/(mol Ang) Convergence criteria 3: must have maximum cell parameter gradient less than this value for convergence
Gamma_Init = 1; % Initial step size coefficient (multiplied by numerical derivative along each dimension) used for moving along gradient
Gamma_Multiplier = 2; %1.5; % Scale Gamma by this value after each cycle (should be greater than 1)
Max_Step_size = 0.5; % Maximum possible step size in angstroms
alpha = 1/2;
beta_JC = 0.5; % Reduce Gamma by this multiple when step size is too large (for JC models).
beta_TF = 0.1; % Reduce Gamma by this multiple when step size is too large (for TF models).
MaxLineTries = 10; % Maximum number of tries to compute lower energy point before decreasing numerical derivative step size
Max_cycle_restarts = 10; % Maximum number of cycle restarts
OptPos_SingleStage = false; % [Not implemented yet] When true, optimize positions and lattice parameters simultaneously by numerical gradients
E_Unphys = -2.2e3; % unphysical energy cutoff
DelE_Unphys = 200; % Unphysical gradient cutoff

%% TOPOLOGY GLOBAL SETTINGS
Top.gen_pairs = 'no'; % Automatically generate pairs
Top.fudgeLJ = 1.0; % Rescale LJ interaction by this amount for 1-4 bonded atoms
Top.fudgeQQ = 1.0; % Rescale Coulomb interaction by this amount for 1-4 bonded atoms

%% MDP GLOBAL SETTINGS
MDP.nsteps_point = 0; % Number of steps to perform before ending (should be 0 for single point energy calculations)
MDP.point_integrator = 'md'; % What type of calculation is run for single point energy calculations (steep = energy min, md = molecular dynamics)
MDP.dt = 0.002; % Time step in ps for md type calculations
MDP.min_integrator = 'steep'; % 'steep' or 'l-bfgs'
MDP.nsteps_min = 1000; % Number of steps to perform before stopping for energy minimization runs
MDP.LJtol = 1e-5; % When doing PME for VdW-interactions, this is used to control the relative strength of the dispersion potential at rvdw in the same way as ewald-rtol controls the electrostatic potential.
MDP.CutOffScheme = 'Verlet'; % Either 'group' or 'Verlet' (does NOT apply to tabulated potentials, these are set to group)
MDP.VerletBT = -1; % This sets the maximum allowed error for pair interactions per particle caused by the Verlet buffer, which indirectly sets rlist unless set to -1, in which case rlist will be used.
MDP.VDWType = 'Cut-off'; % Define the type of van der waals potential used. One of 'PME', 'Cut-off', 'Shift', or 'Switch'
MDP.Auto_Cutoff = false; % Set auto or manual cutoff distances in nm (true = automatically half the smallest box length)
MDP.RList_Cutoff = 2.0; % nm. This should be larger or equal to RCoulomb/RVDW
MDP.RCoulomb_Cutoff = 2.0; % nm. if set to less than 0, then Rc = a;
MDP.RVDW_Cutoff = 2.0; % nm. note that rlist ? rCoulomb = RVDW when using Verlet and VerletBT = -1
MDP.rvdw_switch = 1.9; % nm. where to start switching the LJ potential. Only applies when VDWType = switch
MDP.Fourier_Spacing = 0.10; % nm. Default 0.12 nm. Grid dimensions in PME are controlled with fourierspacing
MDP.PME_Order = 4; % Interpolation order for PME. 4 equals cubic interpolation (default).
MDP.Ewald_rtol = 1e-7; %1e-7 Default (1e-5) The relative strength of the Ewald-shifted direct potential at rcoulomb. Decreasing this will give a more accurate direct sum, but then you need more wave vectors for the reciprocal sum.
MDP.vdw_modifier = 'Potential-shift'; % Potential-shift-Verlet, Potential-shift, None, Force-switch, Potential-switch

%% GEOMETRY AND CALCULATION SETTINGS
% Choose initial conditions if no minimum energy exists yet
% Lengths in Angstroms
Cry.BetaBeO.a = 6.1150; % 5.8
Cry.BetaBeO.c = 3.4351; %(sqrt(3)/2).*Cry.BetaBeO.a;
Cry.BetaBeO.x = 0.336; %0.336; % Experimental value
Cry.BetaBeO.y = 0.310; %0.310; % Experimental value
Cry.BetaBeO.FC_Metal  = [Cry.BetaBeO.x 1-Cry.BetaBeO.x 0.000];
Cry.BetaBeO.FC_Halide = [Cry.BetaBeO.y Cry.BetaBeO.y   0.000];
Cry.BetaBeO.b = Cry.BetaBeO.a;
Cry.BetaBeO.alpha = 90;
Cry.BetaBeO.beta = 90;
Cry.BetaBeO.gamma = 90;
Cry.BetaBeO.Transform = eye(3);
[Cry.BetaBeO.FC_Metal,Cry.BetaBeO.FC_Halide] = ...
    UnitCell_FractionalCoords(Cry.BetaBeO.FC_Metal,...
    Cry.BetaBeO.FC_Halide,'BetaBeO');
Cry.BetaBeO.N = 8;
Cry.BetaBeO.Label = 'B';

Cry.CsCl.a = 3.0343;
Cry.CsCl.FC_Metal  = [0.0 0.0 0.0];
Cry.CsCl.FC_Halide = [1/2 1/2 1/2];
Cry.CsCl.b = Cry.CsCl.a;
Cry.CsCl.c = Cry.CsCl.a;
Cry.CsCl.alpha = 90;
Cry.CsCl.beta = 90;
Cry.CsCl.gamma = 90;
Cry.CsCl.Transform = eye(3);
[Cry.CsCl.FC_Metal,Cry.CsCl.FC_Halide] = ...
    UnitCell_FractionalCoords(Cry.CsCl.FC_Metal,...
    Cry.CsCl.FC_Halide,'CsCl');
Cry.CsCl.N = 2;
Cry.CsCl.Label = 'C';

Cry.FiveFive.a = 4.7273;
Cry.FiveFive.b = 6.7674; %(3/2).*Cry.FiveFive.a; % Theoretical b/a = 3/2
Cry.FiveFive.c = 3.8680; %(sqrt(3)/2).*Cry.FiveFive.a; % Theoretical c/a = cosd(30)
Cry.FiveFive.FC_Metal  = [1/4 1/6 1/2];
Cry.FiveFive.FC_Halide = [1/4 1/3 0.0];
Cry.FiveFive.alpha = 90;
Cry.FiveFive.beta = 90;
Cry.FiveFive.gamma = 90;
Cry.FiveFive.Transform = eye(3);
[Cry.FiveFive.FC_Metal,Cry.FiveFive.FC_Halide] = ...
    UnitCell_FractionalCoords(Cry.FiveFive.FC_Metal,...
    Cry.FiveFive.FC_Halide,'FiveFive');
Cry.FiveFive.N = 8;
Cry.FiveFive.Label = 'F';

Cry.NiAs.a = 2.5339;
Cry.NiAs.x = 1.391666666667;
Cry.NiAs.c = 5.6719; %Cry.NiAs.x.*Cry.NiAs.a; %1.7316 LiF LSDA, 1.7275 LiF PBE
Cry.NiAs.FC_Metal  = [0.0 0.0 0.0];
Cry.NiAs.FC_Halide = [1/3 2/3 1/4];
Cry.NiAs.b = Cry.NiAs.a;
Cry.NiAs.alpha = 90;
Cry.NiAs.beta = 90;
Cry.NiAs.gamma = 120;
Cry.NiAs.Transform =  [1        0        0; ...
                      -sind(30) cosd(30) 0; ...
                       0        0        1];
[Cry.NiAs.FC_Metal,Cry.NiAs.FC_Halide] = ...
    UnitCell_FractionalCoords(Cry.NiAs.FC_Metal,...
    Cry.NiAs.FC_Halide,'NiAs');
Cry.NiAs.N = 4;
Cry.NiAs.Label = 'N';

Cry.Rocksalt.a = 4.6170;
Cry.Rocksalt.FC_Metal  = [0.0 0.0 0.0];
Cry.Rocksalt.FC_Halide = [1/2 1/2 1/2];
Cry.Rocksalt.b = Cry.Rocksalt.a;
Cry.Rocksalt.c = Cry.Rocksalt.a;
Cry.Rocksalt.alpha = 90;
Cry.Rocksalt.beta = 90;
Cry.Rocksalt.gamma = 90;
Cry.Rocksalt.Transform = eye(3);
[Cry.Rocksalt.FC_Metal,Cry.Rocksalt.FC_Halide] = ...
    UnitCell_FractionalCoords(Cry.Rocksalt.FC_Metal,...
    Cry.Rocksalt.FC_Halide,'Rocksalt');
Cry.Rocksalt.N = 8;
Cry.Rocksalt.Label = 'R';

Cry.Sphalerite.a = 4.8764;
Cry.Sphalerite.FC_Metal  = [0.0 0.0 0.0];
Cry.Sphalerite.FC_Halide = [1/4 1/4 1/4];
Cry.Sphalerite.b = Cry.Sphalerite.a;
Cry.Sphalerite.c = Cry.Sphalerite.a;
Cry.Sphalerite.alpha = 90;
Cry.Sphalerite.beta = 90;
Cry.Sphalerite.gamma = 90;
Cry.Sphalerite.Transform = eye(3);
[Cry.Sphalerite.FC_Metal,Cry.Sphalerite.FC_Halide] = ...
    UnitCell_FractionalCoords(Cry.Sphalerite.FC_Metal,...
    Cry.Sphalerite.FC_Halide,'Sphalerite');
Cry.Sphalerite.N = 8;
Cry.Sphalerite.Label = 'S';


Cry.Wurtzite.a = 3.4732;
Cry.Wurtzite.c = 5.2644; %sqrt(8/3).*Cry.Wurtzite.a; % Perfect Wurtzite c/a = sqrt(8/3);
Cry.Wurtzite.FC_Metal  = [1/3 2/3 0.0]; %[1/3 2/3 0.0];
Cry.Wurtzite.FC_Halide = [1/3 2/3 0.6048]; %[1/3 2/3 3/8];(0.3333, 0.6667, 0.9899) (0.3333, 0.6667, 0.3851)
Cry.Wurtzite.b = Cry.Wurtzite.a;
Cry.Wurtzite.alpha = 90;
Cry.Wurtzite.beta = 90;
Cry.Wurtzite.gamma = 120;
Cry.Wurtzite.Transform =  [1        0        0; ...
                          -sind(30) cosd(30) 0; ...
                           0        0        1];
[Cry.Wurtzite.FC_Metal,Cry.Wurtzite.FC_Halide] = ...
    UnitCell_FractionalCoords(Cry.Wurtzite.FC_Metal,...
    Cry.Wurtzite.FC_Halide,'Wurtzite');
Cry.Wurtzite.N = 4;
Cry.Wurtzite.Label = 'W';

%% ADDITIONAL SETTINGS
Settings.Hours = 24; % Max time for job (hours)
Settings.Mins = 0; % Max time for job (minutes)
Settings.nMols_per_Task = -1; % -1 to fix the number of cores
Settings.nCores = 1; % Number of cores to request for calculation (currently limited to 1)
Settings.nTasks_per_Node = 1; % Cores per node to request
Settings.Mempernode = '-1'; % Memory request for server (default = '-1', max per core = '0', use '3G' for Cedar)
Settings.Table_Length = max([MDP.RList_Cutoff MDP.RCoulomb_Cutoff MDP.RVDW_Cutoff])+1.01; % How far should non-automatic tables extend in nm
Settings.Delete_MDPout = true; % Automatically delete MDP out files if true
Settings.Delete_MDlog = true; % Delete log files if true
Settings.Delete_ConfOut = false; % Delete the output configuration if true
Settings.Delete_TRR = false; % Delete the trr file after used if true
Settings.Delete_TPR = true; % Delete the tpr file after used if true
Settings.Delete_EnergyIn = true; % Delete the binary energy file after used if true
Settings.Delete_Supercell = false; % Automatically delete the supercell file if true
Settings.Delete_Backups = true; % Automatically delete the gromacs backup files if true
Settings.CoordType = 'g96'; % Either pdb, gro, or g96 (extra precision)
Settings.Tab_StepSize = 0.0005; % Step size of tabulated potentials in nm

%% Error checks
% Parameter matrix input check
[Row_Par,Col_Par] = size(Parameters);

% Model types
if ~ischar(Model)
    error('Input variable ''Model'' must be a character array')
end

if strcmp(Model,'JC')
    if Row_Par ~= 2 || (Col_Par ~= 2 && Col_Par ~= 3)
        error(['JC model parameter matrix must be 2 x 2 or 2 x 3. ' ...
            'Current input matrix is ' num2str(Row_Par) ' x ' num2str(Col_Par)]);
    end
elseif strcmp(Model,'TF')
    
    if Row_Par ~= 4 || Col_Par ~= 3
        error(['TF model parameter matrix must be 4 x 3. ' ...
            'Current input matrix is ' num2str(Row_Par) ' x ' num2str(Col_Par)]);
    end
else
    error(['Unknown Input Model: ' Model])
end

% Structure
if ~ischar(Structure)
    error('Input variable ''Structure'' must be a character array')
end

% Salt/ion types
if ~ischar(Salt)
    error('Input variable ''Salt'' must be a character array')
end

%% Attempt to build a matlab-gromacs interface
Longest_Cutoff = max([MDP.RList_Cutoff MDP.RCoulomb_Cutoff MDP.RVDW_Cutoff]);
home = fileparts(mfilename('fullpath'));
if ispc % for testing
    Tempdir = tempname;
    gmx = 'wsl source ~/.bashrc; gmx_d';
    threadlock = ' -nt 1 -ntmpi 1 -ntomp 1';
elseif isunix
    [~,Servertxt] = system('hostname -s | cut -c 1-3');
    Server = strtrim(Servertxt);
    if strcmp(Server,'ced') || strcmp(Server,'cdr')
        scrdir = getenv('SCRATCH');
        Tempdir = tempname(scrdir);
        gmx = 'gmx_d';
        threadlock = '';
        
        % Error check: is gmx_d loaded? If no, try to load gmx
        [~,gmxtest] = system('command -v gmx_d');
        if isempty(gmxtest)
            [~,gmxtest] = system('command -v gmx');
            if isempty(gmxtest)
                error('neither gmx_d nor gmx are not currently loaded')
            else
                warning('gmx_d not loaded, reverting to gromacs single precision')
                gmx = 'gmx';
            end
        end
        
    elseif strcmp(Server,'pat')
        Tempdir = tempname;
        gmx = 'source /home/user/Documents/MATLAB/.matlabrc; gmx_d';
        threadlock = ' -nt 1 -ntmpi 1 -ntomp 1';
    elseif strcmpi(Server,'Han') || strcmp(Server,'dhc')
        Tempdir = tempname;
        gmx = 'source ~/.matlabrc; gmx_d';
        threadlock = ' -nt 1 -ntmpi 1 -ntomp 1';
    end
else
    try
        scrdir = getenv('SCRATCH');
        Tempdir = tempname(scrdir);
        gmx = 'gmx_d';
        threadlock = '';
        % Error check: is gmx_d loaded? If no, try to load gmx
        [~,gmxtest] = system('command -v gmx_d');
        if isempty(gmxtest)
            [~,gmxtest] = system('command -v gmx');
            if isempty(gmxtest)
                error('neither gmx_d nor gmx are not currently loaded')
            else
                warning('gmx_d not loaded, reverting to gromacs single precision')
                gmx = 'gmx';
            end
        end
    catch
        error('Unknown Server, manually reset sever settings');
    end
end

if Parallel_Mode
    setenv('OMP_NUM_THREADS','1'); %#ok<*UNRCH>
    setenv('GMX_PME_NUM_THREADS','1');
    setenv('GMX_PME_NTHREADS','1');
    setenv('GMX_OPENMP_MAX_THREADS','1');
    setenv('KMP_AFFINITY','disabled');
    pin = [' -pin on' threadlock];
else
    pin = '';
end

% Load topology template location
Topology_Template = fullfile(home,'templates','Gromacs_Templates',...
'Topology.template');

% Load Topology template
Topology_Template = fileread(Topology_Template);

% Add in global parameters to Topology template
Topology_Template = strrep(Topology_Template,'##GENPAIRS##',Top.gen_pairs);
Topology_Template = strrep(Topology_Template,'##FUDGELJ##',num2str(Top.fudgeLJ));
Topology_Template = strrep(Topology_Template,'##FUDGEQQ##',num2str(Top.fudgeQQ));

% Load mdp template location
MDP.Template = fullfile(home,'templates','Gromacs_Templates',...
'MDP.template');

% Load MDP template
MDP.Template = fileread(MDP.Template);

% Add in global parameters to MDP template
MDP.Template = strrep(MDP.Template,'##NSTEPS##',pad(num2str(MDP.nsteps_point),18));
MDP.Template = strrep(MDP.Template,'##INTEGR##',pad(MDP.point_integrator,18));
MDP.Template = strrep(MDP.Template,'##TIMEST##',pad(num2str(MDP.dt),18));

% Determine type of optimization
if OptPos
    OptTxt = 'FULLOPT';
else
    OptTxt = 'CELLOPT';
end

if ~Maintain_Symmetry
    OptTxt = [OptTxt '_SG1'];
end

% Get Metal and Halide info from Current Salt
[Metal,Halide] = Separate_Metal_Halide(Salt);
Metal_Info = elements('Sym',Metal);
Halide_Info = elements('Sym',Halide);

% Update Directory
Current_Directory = fullfile(Tempdir,Salt);

% Create directory if it does not exist
if ~exist(Current_Directory,'dir')
    mkdir(Current_Directory)
end

% Copy Topology Template
Topology_Temp_Ions = Topology_Template;

% Insert element info into topology template
Topology_Temp_Ions = strrep(Topology_Temp_Ions,'##MET##',pad(Metal,2));
Topology_Temp_Ions = strrep(Topology_Temp_Ions,'##METZ##',pad(num2str(Metal_Info.atomic_number),3));
Topology_Temp_Ions = strrep(Topology_Temp_Ions,'##METMASS##',pad(num2str(Metal_Info.atomic_mass),7));

Topology_Temp_Ions = strrep(Topology_Temp_Ions,'##HAL##',pad(Halide,2));
Topology_Temp_Ions = strrep(Topology_Temp_Ions,'##HALZ##',pad(num2str(Halide_Info.atomic_number),3));
Topology_Temp_Ions = strrep(Topology_Temp_Ions,'##HALMASS##',pad(num2str(Halide_Info.atomic_mass),7));

% Copy MDP Template
MDP.Temp_Ions = MDP.Template;

% Insert salt components into MDP template
MDP.Temp_Ions = strrep(MDP.Temp_Ions,'##MET##',Metal);
MDP.Temp_Ions = strrep(MDP.Temp_Ions,'##HAL##',Halide);

% Copy Topology template
Topology_Temp_Struct = Topology_Temp_Ions;

% Get current structure
Label = Cry.(Structure).Label;

% Select symmetry settings
if Maintain_Symmetry
    switch Structure
        case 'BetaBeO'
            DOF = {'a' 'c'};
            DOF_Units = {[' ' char(0197)] [' ' char(0197)]};
        case 'CsCl'
            DOF = {'a'};
            DOF_Units = {[' ' char(0197)]};
        case 'FiveFive'
            DOF = {'a' 'b' 'c'};
            DOF_Units = {[' ' char(0197)] [' ' char(0197)] [' ' char(0197)]};
        case 'Sphalerite'
            DOF = {'a'};
            DOF_Units = {[' ' char(0197)]};
        case 'NiAs'
            DOF = {'a' 'c'};
            DOF_Units = {[' ' char(0197)] [' ' char(0197)]};
        case 'Rocksalt'
            DOF = {'a'};
            DOF_Units = {[' ' char(0197)]};
        case 'Wurtzite'
            DOF = {'a' 'c'};
            DOF_Units = {[' ' char(0197)] [' ' char(0197)]};
    end
    DOF_txt = DOF;
else
    DOF = {'a' 'b' 'c' 'alpha' 'beta' 'gamma'};
    DOF_txt = {'a' 'b' 'c' char(945) char(946) char(947)};
    DOF_Units = {[' ' char(0197)] [' ' char(0197)] [' ' char(0197)] ...
        [' ' char(0176)] [' ' char(0176)] [' ' char(0176)]};
end

% How many Degrees of freedom
N_DOF = length(DOF);

% Update Directory
Current_Directory = fullfile(Tempdir,Salt,Structure);

% Create directory if it does not exist
if ~exist(Current_Directory,'dir')
    mkdir(Current_Directory)
end

% Structure Template filename
Coordinate_File = fullfile(home,'templates',...
    [upper(Settings.CoordType) '_Templates'],...
    [Structure '.' Settings.CoordType]);

% Load Template text
Coordinate_text = fileread(Coordinate_File);

% Add Metal and Halide symbols
if strcmp(Settings.CoordType,'gro')
    Met = pad(Metal,2,'left');
    Hal = pad(Halide,2,'left');
elseif strcmp(Settings.CoordType,'g96')
    Met = pad(Metal,2,'right');
    Hal = pad(Halide,2,'right');
end
Coordinate_text = strrep(Coordinate_text,'##MET##',Met);
Coordinate_text = strrep(Coordinate_text,'##HAL##',Hal);

% Calculate size of supercell
N_Supercell = ceil((N_atoms/Cry.(Structure).N)^(1/3));

% Add number of unit cells to topology file
Topology_Temp_Struct = strrep(Topology_Temp_Struct,'##N##',num2str(N_Supercell));

% Update Directory
Current_Directory = fullfile(Tempdir,Salt,...
    Structure,Model);

% Create directory if it does not exist
if ~exist(Current_Directory,'dir')
    mkdir(Current_Directory)
end

% Load topology template and MDP template
Topology_text = Topology_Temp_Struct;
MDP.Temp_Model = MDP.Temp_Ions;

% Update Topology and MDP files
MDP.Temp_Model = strrep(MDP.Temp_Model,'##COULOMB##',pad('PME',18));
MDP.Temp_Model = strrep(MDP.Temp_Model,'##FOURIER##',pad(num2str(MDP.Fourier_Spacing),18));
MDP.Temp_Model = strrep(MDP.Temp_Model,'##PMEORDER##',pad(num2str(MDP.PME_Order),18));
MDP.Temp_Model = strrep(MDP.Temp_Model,'##EWALDTOL##',pad(num2str(MDP.Ewald_rtol),18));

% TF model defined as a table in GROMACS
if strcmp(Model,'TF')

    beta = beta_TF;

    % Define the function type as 1 (needed for custom functions)
    Topology_text = strrep(Topology_text,'##NBFUNC##','1');

    % Define the combination rules (Lorenz-berthelot)
    Topology_text = strrep(Topology_text,'##COMBR##','1');

    % Define all the parameters as 1.0 (already included in potentials)
    Topology_text = strrep(Topology_text,'##METMETC##',pad('1.0',10));
    Topology_text = strrep(Topology_text,'##HALHALC##',pad('1.0',10));
    Topology_text = strrep(Topology_text,'##METHALC##',pad('1.0',10));
    Topology_text = strrep(Topology_text,'##METMETA##','1.0');
    Topology_text = strrep(Topology_text,'##HALHALA##','1.0');
    Topology_text = strrep(Topology_text,'##METHALA##','1.0');

    % Generate tables of the TF potential
    [TF_U_PM, TF_U_PP, TF_U_MM] = TF_Potential_Generator(0,...
        Settings.Table_Length,Settings.Tab_StepSize,Salt,Parameters,false,...
        MDP.vdw_modifier,MDP.RVDW_Cutoff);

    TableName = [Salt '_' Model '_Table'];
    TableFile = fullfile(Current_Directory,[TableName '.xvg']);

    % Save tables into current directory
    fidPM = fopen(fullfile(Current_Directory,[TableName '.xvg']),'wt');
    fwrite(fidPM,regexprep(TF_U_PM,{'\r', '\n\n+'}',{'', '\n'}));
    fclose(fidPM);

    fidPP = fopen(fullfile(Current_Directory,...
        [TableName '_' Metal '_' Metal '.xvg']),'wt');
    fwrite(fidPP,regexprep(TF_U_PP,{'\r', '\n\n+'}',{'', '\n'}));
    fclose(fidPP);

    fidMM = fopen(fullfile(Current_Directory,...
        [TableName '_' Halide '_' Halide '.xvg']),'wt');
    fwrite(fidMM,regexprep(TF_U_MM,{'\r', '\n\n+'}',{'', '\n'}));
    fclose(fidMM);

    % Modify the MDP file
    MDP.Temp_Model = strrep(MDP.Temp_Model,'##VDWTYPE##',pad('user',18));
    MDP.Temp_Model = strrep(MDP.Temp_Model,'##VDWMOD##',pad(MDP.vdw_modifier,18));
    MDP.Temp_Model = strrep(MDP.Temp_Model,'##CUTOFF##',pad('group',18));
    MDP.Temp_Model = regexprep(MDP.Temp_Model,'ewald-rtol-lj.+?\n','');
    MDP.Temp_Model = regexprep(MDP.Temp_Model,'lj-pme-comb-rule.+?\n','');
    MDP.Temp_Model = regexprep(MDP.Temp_Model,'verlet-buffer-tolerance.+?\n','');

    % Energy conversion setting
    EnergySetting = '1 2 3 4 28 29 30 31 32 33 0';

    % JC model defined as parameters in GROMACS (faster)
elseif strcmp(Model,'JC') && ~any(Parameters(:) <= 0) && Col_Par == 2

    beta = beta_JC;

    TableFile = '';

    % Definte the function type as 1 (LJ)
    Topology_text = strrep(Topology_text,'##NBFUNC##','1');

    % Define the combination rules (Lorenz-berthelot in sigma-epsilon form)
    Topology_text = strrep(Topology_text,'##COMBR##','2');

    % Get JC parameters
    Met_JC_Param.sigma = Parameters(1,1);
    Hal_JC_Param.sigma = Parameters(1,2);
    Met_JC_Param.epsilon = Parameters(2,1);
    Hal_JC_Param.epsilon = Parameters(2,2);

    % Cross terms using Lorenz-Berthelot combining rules
    Sigma_ij = (1/2)*(Met_JC_Param.sigma + Hal_JC_Param.sigma);
    epsilon_ij = sqrt(Met_JC_Param.epsilon*Hal_JC_Param.epsilon);

    % Add parameters to topology text
    Topology_text = strrep(Topology_text,'##METMETC##',pad(num2str(Met_JC_Param.sigma,'%10.8f'),10));
    Topology_text = strrep(Topology_text,'##HALHALC##',pad(num2str(Hal_JC_Param.sigma,'%10.8f'),10));
    Topology_text = strrep(Topology_text,'##METHALC##',pad(num2str(Sigma_ij,'%10.8f'),10));
    Topology_text = strrep(Topology_text,'##METMETA##',num2str(Met_JC_Param.epsilon,'%10.8f'));
    Topology_text = strrep(Topology_text,'##HALHALA##',num2str(Hal_JC_Param.epsilon,'%10.8f'));
    Topology_text = strrep(Topology_text,'##METHALA##',num2str(epsilon_ij,'%10.8f'));

    % Modify the MDP file
    MDP.Temp_Model = strrep(MDP.Temp_Model,'##VDWTYPE##',pad(MDP.VDWType,18));
    MDP.Temp_Model = strrep(MDP.Temp_Model,'##VDWMOD##',pad(MDP.vdw_modifier,18));
    MDP.Temp_Model = strrep(MDP.Temp_Model,'##CUTOFF##',pad(MDP.CutOffScheme,18));
    MDP.Temp_Model = regexprep(MDP.Temp_Model,'energygrp-table.+?\n','');
    MDP.Temp_Model = regexprep(MDP.Temp_Model,'ewald-rtol-lj.+?\n','');
    MDP.Temp_Model = regexprep(MDP.Temp_Model,'lj-pme-comb-rule.+?\n','');

    % Add in Verlet Settings
    if strcmp(MDP.CutOffScheme,'Verlet')
        MDP.Temp_Model = strrep(MDP.Temp_Model,'##VerletBT##',pad(num2str(MDP.VerletBT),18));
    else
        MDP.Temp_Model = regexprep(MDP.Temp_Model,'verlet-buffer-tolerance.+?\n','');
    end

    % Energy conversion setting
    EnergySetting = '1 2 3 4 28 29 30 31 32 33 0';

    % JC model defined as table in GROMACS (slower)
elseif strcmp(Model,'JC') && (any(Parameters <= 0,'all') || Col_Par == 3)

    beta = beta_JC;

    % Define the function type as 1 (needed for custom functions)
    Topology_text = strrep(Topology_text,'##NBFUNC##','1');

    % Define the combination rules (not used for tables input)
    Topology_text = strrep(Topology_text,'##COMBR##','1');

    % Define all the parameters as 1.0 (already included in potentials)
    Topology_text = strrep(Topology_text,'##METMETC##',pad('1.0',10));
    Topology_text = strrep(Topology_text,'##HALHALC##',pad('1.0',10));
    Topology_text = strrep(Topology_text,'##METHALC##',pad('1.0',10));
    Topology_text = strrep(Topology_text,'##METMETA##','1.0');
    Topology_text = strrep(Topology_text,'##HALHALA##','1.0');
    Topology_text = strrep(Topology_text,'##METHALA##','1.0');

    % Generate tables of the JC potential
    [JC_U_PM, JC_U_PP, JC_U_MM] = JC_Potential_Generator(0,...
        Settings.Table_Length,Settings.Tab_StepSize,Salt,Parameters,false,...
        MDP.vdw_modifier,MDP.RVDW_Cutoff);            

    TableName = [Salt '_' Model '_Table'];
    TableFile = fullfile(Current_Directory,[TableName '.xvg']);

    % Save tables into current directory
    fidPM = fopen(fullfile(Current_Directory,[TableName '.xvg']),'wt');
    fwrite(fidPM,regexprep(JC_U_PM,{'\r', '\n\n+'}',{'', '\n'}));
    fclose(fidPM);

    fidPP = fopen(fullfile(Current_Directory,...
        [TableName '_' Metal '_' Metal '.xvg']),'wt');
    fwrite(fidPP,regexprep(JC_U_PP,{'\r', '\n\n+'}',{'', '\n'}));
    fclose(fidPP);

    fidMM = fopen(fullfile(Current_Directory,...
        [TableName '_' Halide '_' Halide '.xvg']),'wt');
    fwrite(fidMM,regexprep(JC_U_MM,{'\r', '\n\n+'}',{'', '\n'}));
    fclose(fidMM);

    % Modify the MDP file
    MDP.Temp_Model = strrep(MDP.Temp_Model,'##VDWTYPE##',pad('user',18));
    MDP.Temp_Model = strrep(MDP.Temp_Model,'##VDWMOD##',pad(MDP.vdw_modifier,18));
    MDP.Temp_Model = strrep(MDP.Temp_Model,'##CUTOFF##',pad('group',18));
    MDP.Temp_Model = regexprep(MDP.Temp_Model,'ewald-rtol-lj.+?\n','');
    MDP.Temp_Model = regexprep(MDP.Temp_Model,'lj-pme-comb-rule.+?\n','');
    MDP.Temp_Model = regexprep(MDP.Temp_Model,'verlet-buffer-tolerance.+?\n','');

    % Energy conversion setting
    EnergySetting = '1 2 3 4 28 29 30 31 32 33 0';
end

TotalTimer = tic;
disp(['Beginning ' Salt ' ' Structure ' ' Model ' Optimization...'])
disp(['Model Chosen: ' Model])
disp('Parameter input: ')
for i = 1:size(Parameters,1)
    disp(num2str(Parameters(i,:),'%1.4E '))
end

% Update directory details
FileBase = [Salt '_' Label '_' Model '_' OptTxt];

if ispc
    OptDir = [Current_Directory filesep OptTxt];
    OptDirUnix = windows2unix([Current_Directory filesep OptTxt]);
    rm_command = ['wsl find ' OptDirUnix ' -iname "#*#" ^| xargs rm -f'];
    del_command = ['wsl rm -r ' OptDirUnix '/*'];
else
    OptDir = [Current_Directory filesep OptTxt];
    rm_command = ['find ' OptDir ' -iname "#*#" | xargs rm -f'];
    del_command = ['rm -r ' OptDir '/*'];
end

% Find minimum lattice parameter for this
% salt/structure/model (or use initial ones)
if Find_Min_Params
    [Cry.(Structure).a,Cry.(Structure).b,...
        Cry.(Structure).c,...
        Cry.(Structure).FC_Metal,...
        Cry.(Structure).FC_Halide,Updated] = ...
        FindMinLatticeParam(Cry,Salt,...
        Structure,Model,home,Data_Types,...
        Find_Similar_Params);

    if ~Continue_if_no_IC && ~Updated
        disp(['No Suitable Initial Conditions Found. For '  Salt ' ' Structure ' ' Model ' Optimization.'])
        disp('Calculation Halted.')
        system(del_command);
         % Delete tables
        if strcmp(Model,'TF') || (strcmp(Model,'JC') && Dispersion_Scale <= 0)
            delete(fullfile(Current_Directory,[TableName '.xvg']));
            delete(fullfile(Current_Directory,...
                [TableName '_' Metal '_' Metal '.xvg']));
            delete(fullfile(Current_Directory,...
                [TableName '_' Halide '_' Halide '.xvg']));
        end
        return
    end
end

%% Begin Optimization Loop
% Loop broken when convergence criteria is met or max cycles reached
skip_results = false;
auto_h = true;
[~,~] = system(rm_command);
Cycle_restarts = 0;
Gamma = Gamma_Init;
for Index = 1:MaxCycles
    Restart_cycle = false;
    % Calculate reasonable step size
    if auto_h
        h = nuderst([Cry.(Structure).a Cry.(Structure).b Cry.(Structure).c ...
            Cry.(Structure).alpha Cry.(Structure).beta Cry.(Structure).gamma]);
        
        OptimizationLoop(gmx,Cry,Salt,Structure,Model,Label,N_Supercell,...
            Tempdir,MDP,Longest_Cutoff,Coordinate_text,Settings,...
            Topology_text,TableFile,EnergySetting,OptTxt,pin);
        E = GrabEnergy(OptDir,FileBase);

        disp(['********************Cycle ' num2str(Index) ' Initial Conditions********************']);
        for Didx = 1:N_DOF
            disp(['Lattice Parameter ' DOF_txt{Didx} ' = ' ...
                num2str(Cry.(Structure).(DOF{Didx}),'%2.8f') DOF_Units{Didx} '.' ...
                ' Num. Der. Step Size: ' char(948) '(' DOF_txt{Didx} ') = ' ...
                num2str(h(Didx),'%2.8f') DOF_Units{Didx} '.']);
        end
        disp(['Fractional Coordinates for ' Metal ': ']);
        disp(num2str(Cry.(Structure).FC_Metal(:,:),'%2.8f '))
        disp(['Fractional Coordinates for ' Halide ': ']);
        disp(num2str(Cry.(Structure).FC_Halide(:,:),'%2.8f '))
        disp(['Initial E = ' num2str(E,'%4.10f') ' kJ/mol']);
        disp(['Step size coefficient: ' num2str(Gamma) ])
        Telap = datestr(seconds(toc(TotalTimer)),'HH:MM:SS');
        disp(['Time Elapsed: ' Telap])
        disp('******************************************************************')
    else
        auto_h = true;
    end

    % Check for unphysical energy
    if E < E_Unphys
        disp(['Warning: Unphysical energy detected after ' num2str(Index) ' cycles.'])
        disp('Model may have local minimum at complete overlap of opposite charges.')
        disp('Removing output files.')
        system(del_command);
        skip_results = true;
        E = nan;
        Cry.(Structure).a = nan;
        Cry.(Structure).b = nan;
        Cry.(Structure).c = nan;
        Cry.(Structure).FC_Metal(:) = nan;
        Cry.(Structure).FC_Halide(:) = nan;
        break
    end
      
    % Loop through each DOF
    Gradient = nan(1,N_DOF);
    for Didx = 1:N_DOF

        CurDOF = DOF{Didx};
        CurDOFtxt = DOF_txt{Didx};
        CurDOFunit = DOF_Units{Didx};
        delta = h(Didx);

        % Get energy at minus step
        CryMinus = Cry;
        CryMinus.(Structure).(CurDOF) = CryMinus.(Structure).(CurDOF) - delta;

        % Update transformation matrix
        CryMinus.(Structure).Transform = GenTransformMatrix(CryMinus.(Structure));

        % Maintain symmetry if set
        if Maintain_Symmetry && strcmp(CurDOF,'a')
            [CryMinus.(Structure).b,CryMinus.(Structure).c] = SymmetryAdapt(CryMinus.(Structure).a,...
                CryMinus.(Structure).b,CryMinus.(Structure).c,Structure); %#ok<*UNRCH>
        end

        % Generate energy at new point
        OptimizationLoop(gmx,CryMinus,Salt,Structure,Model,Label,N_Supercell,...
            Tempdir,MDP,Longest_Cutoff,Coordinate_text,Settings,...
            Topology_text,TableFile,EnergySetting,OptTxt,pin);
        E_Minus = GrabEnergy(OptDir,FileBase);

        Telap = datestr(seconds(toc(TotalTimer)),'HH:MM:SS');
        disp(['Finite step ' CurDOFtxt ' -' num2str(delta,'%2.2E') ...
            CurDOFunit ' Calculated. ' char(916) 'E(' CurDOFtxt '-' char(948) ...
            ') = ' num2str(E_Minus-E,'%+4.4E') ' kJ/mol. Time Elapsed: ' Telap]);

        % Get energy at plus step
        CryPlus = Cry;
        CryPlus.(Structure).(CurDOF) = CryPlus.(Structure).(CurDOF) + delta;

        % Update transformation matrix
        CryPlus.(Structure).Transform = GenTransformMatrix(CryPlus.(Structure));

        % Maintain symmetry if set
        if Maintain_Symmetry && strcmp(CurDOF,'a')
            [CryPlus.(Structure).b,CryPlus.(Structure).c] = SymmetryAdapt(CryPlus.(Structure).a,...
                CryPlus.(Structure).b,CryPlus.(Structure).c,Structure); %#ok<*UNRCH>
        end

        % Generate Energy at new point
        OptimizationLoop(gmx,CryPlus,Salt,Structure,Model,Label,N_Supercell,...
            Tempdir,MDP,Longest_Cutoff,Coordinate_text,Settings,...
            Topology_text,TableFile,EnergySetting,OptTxt,pin);
        E_Plus = GrabEnergy(OptDir,FileBase);

        Telap = datestr(seconds(toc(TotalTimer)),'HH:MM:SS');
        disp(['Finite step ' CurDOFtxt ' +' num2str(delta,'%2.2E') ...
            CurDOFunit ' Calculated. ' char(916) 'E(' CurDOFtxt '+' char(948) ...
            ') = ' num2str(E_Plus-E,'%+4.4E') ' kJ/mol. Time Elapsed: ' Telap]);

        % 3-Point Numerical derivative wrt current DOF
        Gradient(Didx) = (1/(2*delta))*(-E_Minus + E_Plus);
        [~,~] = system(rm_command);
    end
    disp([char(8711) 'E = ' num2str(Gradient,'%5.4E  ') ' kJ/mol ' char(0197)])
    
    % Check for unphysical negative gradient
    if max(Gradient) > DelE_Unphys
        disp(['Warning: Unphysical energy gradient detected after ' num2str(Index) ' cycles.'])
        disp('Model may have local minimum at complete overlap of opposite charges.')
        disp('Othterwise poor initial conditions may cause this.')
        disp('Removing output files.')
        system(del_command);
        skip_results = true;
        E = nan;
        Cry.(Structure).a = nan;
        Cry.(Structure).b = nan;
        Cry.(Structure).c = nan;
        Cry.(Structure).FC_Metal(:) = nan;
        Cry.(Structure).FC_Halide(:) = nan;
        break
        % Check for convergence on cycle 1
    elseif Index == 1 && (rms(Gradient) < Gradient_Tol_RMS) && ...
            (max(abs(Gradient)) < Gradient_Tol_Max)
        % If gradient convergence criteria are met, end loop
        disp('Energy convergence reached on cycle 1 gradient test.')
        [~,~] = system(rm_command);
        break
    end
    
    % Move one step in direction of steepest descent
    for StepInd = 1:MaxLineTries

        CryNew = Cry;

        for Didx = 1:N_DOF

            CurDOF = DOF{Didx};
            CurGrad = Gradient(Didx);

            CryNew.(Structure).(CurDOF) = Cry.(Structure).(CurDOF)...
                - sign(CurGrad)*min(abs(Gamma*CurGrad),Max_Step_size);
        end
        if Maintain_Symmetry
            [CryNew.(Structure).b,CryNew.(Structure).c] = SymmetryAdapt(CryNew.(Structure).a,...
                CryNew.(Structure).b,CryNew.(Structure).c,Structure);
        end

        % Update transformation matrix
        CryNew.(Structure).Transform = GenTransformMatrix(CryNew.(Structure));

        % Recalculate energy at new point
        OptimizationLoop(gmx,CryNew,Salt,Structure,Model,Label,N_Supercell,...
            Tempdir,MDP,Longest_Cutoff,Coordinate_text,Settings,...
            Topology_text,TableFile,EnergySetting,OptTxt,pin);
        E_New = GrabEnergy(OptDir,FileBase);

        if E_New > E - min(alpha*Gamma*norm(Gradient)^2,50)
            Gamma = beta*Gamma;

            Telap = datestr(seconds(toc(TotalTimer)),'HH:MM:SS');
            disp(['Optimizing step size... Time Elapsed: ' Telap])
            [~,~] = system(rm_command);
        else
            Telap = datestr(seconds(toc(TotalTimer)),'HH:MM:SS');
            disp(['Finite step in direction of steepest descent calculated. Time Elapsed: ' Telap]);
            disp([char(916) 'E = ' num2str(E_New-E,'%4.10f') ' kJ/mol']);
            disp(['RMS ' char(8711) 'E = ' num2str(rms(Gradient)) ' kJ/mol ' char(0197)])
            disp(['|Max(' char(8711) 'E)| = ' num2str(max(abs(Gradient))) ' kJ/mol ' char(0197)])
            break
        end

        if StepInd >= MaxLineTries
            h = max(h./10,eps*100); % dont go below eps
            Gamma = Gamma_Init;
            Restart_cycle = true;
            [~,~] = system(rm_command);
            break
        end
    end
    
    if Restart_cycle
        Cycle_restarts = Cycle_restarts+1;
        if Cycle_restarts >= Max_cycle_restarts
            disp(['Unable to find lower energy point after ' num2str(Max_cycle_restarts) ' attempts.'])
            disp('Stopping here.')
            if (rms(Gradient) < Gradient_Tol_RMS) && (max(abs(Gradient)) < Gradient_Tol_Max)
                disp('Gradient convergence criteria fulfilled.')
                disp(['Energy convergence reached after ' num2str(Index) ' cycles.'])
            else
                disp('Gradient convergence criteria NOT fulfilled.')
                disp(['Search haulted after ' num2str(Index) ' cycles.'])
                disp('Check initial conditions. Poor initial conditions may cause this.')
                disp('Removing output files.')
                system(del_command);
                skip_results = true;
                E = nan;
                Cry.(Structure).a = nan;
                Cry.(Structure).b = nan;
                Cry.(Structure).c = nan;
                Cry.(Structure).FC_Metal(:) = nan;
                Cry.(Structure).FC_Halide(:) = nan;
            end
            break
        else
            auto_h = false;
            disp('Unable to find a point of lower energy in chosen direction, starting next cycle with modified numerical derivative step size.');
            continue
        end
    end
    
    % Check for unphysical energy
    if abs(E_New - E) > 4e3 || (E_New < E_Unphys) || (max(Gradient) > DelE_Unphys)
        disp(['Warning: Unphysical energy detected after ' num2str(Index) ' cycles.'])
        disp('Model may have local minimum at complete overlap of opposite charges.')
        disp('Othterwise poor initial conditions may cause this.')
        disp('Removing output files.')
        system(del_command);
        skip_results = true;
        E = nan;
        Cry.(Structure).a = nan;
        Cry.(Structure).b = nan;
        Cry.(Structure).c = nan;
        Cry.(Structure).FC_Metal(:) = nan;
        Cry.(Structure).FC_Halide(:) = nan;
        break
    end
    
    % Optimize wrt positions
    if OptPos
        % Save MDP variable for later
        Temp_ModelOld = MDP.Temp_Model;

        MDP.Temp_Model = strrep(MDP.Temp_Model,'= md                ; What type of calculation is run',...
            ['= ' MDP.min_integrator newline 'emtol                    = ' num2str(emtol)]);
        MDP.Temp_Model = regexprep(MDP.Temp_Model,'nsteps                   = [0-9|\.|\-]+',...
            ['nsteps                   = ' num2str(MDP.nsteps_min)]);

        % Rerun energy calc and break loop
        EnergySettingAlt = '1 2 3 4 25 26 27 28 29 30 0'; 

        N_Supercell_out = OptimizationLoopFC(gmx,CryNew,Salt,Structure,...
            Model,Label,N_Supercell,...
            Tempdir,MDP,Longest_Cutoff,Coordinate_text,Settings,...
            Topology_text,TableFile,EnergySettingAlt,OptTxt,pin);

        E_Old = E_New;
        E_New = GrabEnergyFinal(OptDir,FileBase);
        disp(['Geometry Optimized W.R.T. Atomic Positions. ' char(916) 'E = ' num2str(E_New-E,'%4.10f') ' kJ/mol']);
        Telap = datestr(seconds(toc(TotalTimer)),'HH:MM:SS');
        disp(['Time Elapsed: ' Telap])

        % Get new unit cell coordinates

        % Get output config filename
        OutConf_file = fullfile(OptDir,[FileBase 'OutConf.g96']);

        % Open the unit cell structure file
        SuperCellText = fileread(OutConf_file);
        SuperCellText = regexprep(SuperCellText, {'\r', '\n\n+'}, {'', '\n'});
        Positions = regexp(SuperCellText,'(?<=POSITION\n)(.+?)(?=(\n)END)','match','ONCE');

        % g96 file format
        Coords_Data = textscan(Positions,'%*6c%6c%*6c%*6c%15.9f%15.9f%15.9f\n',...
            CryNew.(Structure).N,'Delimiter','','whitespace','');

        % Get xyz coordinates of the unit cell
        idx = (CryNew.(Structure).N)/2 + 1;
        Met_XYZ = [Coords_Data{2}(1:idx-1) Coords_Data{3}(1:idx-1) Coords_Data{4}(1:idx-1)];
        Hal_XYZ = [Coords_Data{2}(idx:end) Coords_Data{3}(idx:end) Coords_Data{4}(idx:end)];

        % Get end cell parameter text
        boxtemp = regexp(SuperCellText,'(?<=BOX\n)(.+?)(?=(\n)END)','match');
        boxcoords = textscan(boxtemp{1},'%f %f %f %f %f %f %f %f %f',...
            'Delimiter',' ','MultipleDelimsAsOne',true);

        % Check for empty coords
        indx = cellfun('isempty',boxcoords); % true for empty cells
        boxcoords(indx) = {0}; % replace by a cell with a zero

        % Lattice vectors
        a_vec = [boxcoords{1} boxcoords{4} boxcoords{5}]/N_Supercell_out;
        b_vec = [boxcoords{6} boxcoords{2} boxcoords{7}]/N_Supercell_out;
        c_vec = [boxcoords{8} boxcoords{9} boxcoords{3}]/N_Supercell_out;

        if (strcmp(Structure,'Wurtzite') || strcmp(Structure,'NiAs')) && VecAngle(a_vec,b_vec) < 90.0
            b_vec = [-boxcoords{6} boxcoords{2} boxcoords{7}]/N_Supercell_out;
        end

        % Get transformation Matrix
        TM = [a_vec ; b_vec ; c_vec];
        CryNew.(Structure).Transform = [a_vec/norm(a_vec) ; b_vec/norm(b_vec) ; c_vec/norm(c_vec)];

        % Convert to Fractional coordinates
        CryNew.(Structure).FC_Metal = mod(Met_XYZ/TM,1);
        CryNew.(Structure).FC_Halide = mod(Hal_XYZ/TM,1);

        % Copy to unit cell
        Unit_cell_text = AddCartesianCoord(Coordinate_text,CryNew.(Structure),1,false,'g96');

        % Overwrite unit cell
        UnitCell_file = fullfile(OptDir,[FileBase '_UnitCell.g96']);
        fidUC = fopen(UnitCell_file,'wt');
        fwrite(fidUC,regexprep(Unit_cell_text,{'\r', '\n\n+'}',{'', '\n'}));
        fclose(fidUC);

        % Reload old MDP variable
        MDP.Temp_Model = Temp_ModelOld;
    end

    % Remove backup files
    [~,~] = system(rm_command);

    Telap = datestr(seconds(toc(TotalTimer)),'HH:MM:SS');
    disp('************************************')
    disp(['Cycle ' num2str(Index) ' complete.'])
    disp(['Time Elapsed: ' Telap])
    disp('************************************')
                        
                        
    % Check if energy converged
    if (abs(E_New - E) < Energy_Tol) && ...
            (rms(Gradient) < Gradient_Tol_RMS) && ...
            (max(abs(Gradient)) < Gradient_Tol_Max)
        % If all convergence criteria are met, end loop
        Cry = CryNew;
        disp(['Energy convergence reached after ' num2str(Index) ' cycles.' newline 'Final recalculation of energy...'])
        OptimizationLoop(gmx,Cry,Salt,Structure,Model,Label,N_Supercell,...
            Tempdir,MDP,Longest_Cutoff,Coordinate_text,Settings,...
            Topology_text,TableFile,EnergySetting,OptTxt,pin);

        E = GrabEnergy(OptDir,FileBase);
        [~,~] = system(rm_command);
        break
    % If energy converged but not gradients
    elseif (abs(E_New - E) < Energy_Tol)
        Gamma = Gamma_Init; % Re-initialize Gamma
        E = E_New;
        Cry = CryNew;
    % Detect unphysical energy
    elseif abs(E_New - E) > 4e3 || (E_New < E_Unphys) || (max(Gradient) > DelE_Unphys)
        disp(['Warning: Unphysical energy detected after ' num2str(Index) ' cycles.'])
        disp('Model may have local minimum at complete overlap of opposite charges.')
        disp('Othterwise poor initial conditions may cause this.')
        disp('Removing output files.')
        system(del_command);
        skip_results = true;
        E = nan;
        Cry.(Structure).a = nan;
        Cry.(Structure).b = nan;
        Cry.(Structure).c = nan;
        Cry.(Structure).FC_Metal(:) = nan;
        Cry.(Structure).FC_Halide(:) = nan;
        break
    % Otherwise a normal decrease in energy
    elseif E_New < E
        Gamma = Gamma*Gamma_Multiplier;
        E = E_New;
        Cry = CryNew;
    % If energy increases
    elseif E_New > E
        disp('Warning: Total energy increased after geometry optimization of atomic positions.')
        disp('Convergence may not have been reached by GROMACS atomic position minimization step.')
        disp('Reverting to intial fractional coordinates.');
        Gamma = Gamma*Gamma_Multiplier;
        E = E_Old;
        CryNew.(Structure).FC_Metal = Cry.(Structure).FC_Metal;
        CryNew.(Structure).FC_Halide = Cry.(Structure).FC_Halide;
        Cry = CryNew;
    end

    if Index < MaxCycles
        disp(['Beginning Cycle ' num2str(Index+1) '.'])
    else
        disp(['Convergence NOT reached after ' num2str(Index) ' cycles. Stopping.'])
    end                        
end             

Telap = datestr(seconds(toc(TotalTimer)),'HH:MM:SS');
disp(['Completed: ' Salt ' ' Structure ' ' Model ' Geometry Optimization. Time ' Telap])
if ~skip_results
    disp(['Final Optimized Energy is ' num2str(E,'%4.10f') ' kJ/mol'])
    for Didx = 1:N_DOF
        disp(['Lattice Parameter ' DOF_txt{Didx} ' = ' ...
            num2str(Cry.(Structure).(DOF{Didx}),'%2.8f') DOF_Units{Didx} '.']);
    end
    disp(['Fractional Coordinates for ' Metal ': ']);
    disp(num2str(Cry.(Structure).FC_Metal(:,:),'%2.8f '))
    disp(['Fractional Coordinates for ' Halide ': ']);
    disp(num2str(Cry.(Structure).FC_Halide(:,:),'%2.8f '))
    disp('Components of Final Gradient:')
    disp(num2str(Gradient,'%5.4E  '))
end

%% Package output
Output_Array = [E Cry.(Structure).a Cry.(Structure).b Cry.(Structure).c ...
    Cry.(Structure).FC_Metal(1,:)  Cry.(Structure).FC_Halide(1,:)];

% Cleanup temp directory
try
    rmdir(Tempdir,'s');
catch
    disp(['Warning - unable to delete temporary directory: ' Tempdir])
end
end
