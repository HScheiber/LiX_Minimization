%% Gen_Molten_Salts_Batch
% A script for generating MD inputs for Lithium halides salts using GROMACS

% Damp_Types:
% 0 = no (default) damping. This is default of JC model.
% 1 = BJ/rational damping (same as in D3(BJ))
% 2 = Tang Damping (Essentially removes dispersion in JC)
% 3 = MMDRE Damping function (Very weak damping, damps mainly at mid range)
% 4 = PAMoC Damping function (fairly weak damping, damps mainly at mid range)
% 5 = EHFSK Damping function (strong damping)
% 6 = WY damping function (strongest damping)

% TF Parameter sets for C6/C8 coefficients (for when using as a script)
% 0 = default TF Parameters
% 1 = D3 values
% 2 = Best literature values available
% 3 = D4 with C6 and C8 generated on the fly

% ModelID is a ID to keep track of your input model. Will be added to your
% input model name. Can be a string or a numeric.

% Job name convention: 
%$MAINDIR/Salt/(Date)_(Initial Struc Label)_(Model)_(ModelID)_(NVE or NPT)/(Job Files)

% For inputing your own JC/TF parameters use the following input order:
% Gen_Molten_Salts_Batch(Salt,Structure,Model,Parameters,OptPos,CRDamping,C6Damping,ModelID)

function Gen_Molten_Salts_Batch(varargin)
%#ok<*NASGU>
%#ok<*STRNU>
%% Check inputs
if nargin ~= 0
    % Ensure correct number of inputs
    if nargin == 8
        Salt = varargin{1};
        Structure = varargin{2};
        Model = varargin{3};
        Parameters = varargin{4};
        OptPos = varargin{5};
        CRDamping = varargin{6};
        Damp = varargin{7};
        ModelID = varargin{8};
        Potential_Input = true;
        TF_Param = 0;
    else
        error('Number of inputs should either be 0 or 8. Input order is: (Salt,Structure,Model,Parameters,OptPos,CRDamping,C6Damping,ModelID)');
    end
else
    Potential_Input = false;
end

%% User settings
CC_Username = 'dzhp1'; % Compute Canada Username (used for graham, cedar, orcinus and beluga)
CWL_Username = 'haydensc'; % CWL username (only used for sockeye);
LiX_Directory ='scratch/CHEM449/LiX_Minimization'; % location of directory

%% Job Settings (any of these can be arrays)
if ~Potential_Input
    Salts = {'LiF' 'LiCl' 'LiBr' 'LiI' 'NaCl'};
    Structures = {'Rocksalt'}; %{'Rocksalt' 'Wurtzite' 'Sphalerite' 'CsCl' 'NiAs' 'BetaBeO' 'FiveFive'}; Initial structure
    Models = {'JC'}; % Input model(s) to use: JC, JC3P, JC4P, TF
    Damping_Funcs = 0; % Input damping functions to use: 0 to 6 are defined
    TF_Paramset = 0; % choose from 0-3
    Scale_Dispersion = 1.0; % Works for both JC and TF
    Scale_Repulsion = 1.0; % Works for both JC and TF
    Scale_MM_Dispersion = 1.0; % Works for both JC and TF
    Scale_XX_Dispersion = 1.0; % Works for both JC and TF
    Scale_MX_Dispersion = 1.0; % Works for both JC and TF
    Scale_Epsilon = 1.0; % Scale all Epsilon (affects JC only)
    Scale_Sigma = 1.0; % Scale all Sigma (affects JC only)
    Scale_Alpha = 1.0; % Scale the repulsive exponential parameter alpha (affects TF only)
end

%% Compute Node settings
JobSettings.Hours = 24; % Max time for each job (hours)
JobSettings.Mins = 0; % Max time for job (minutes)
JobSettings.Cores = 32; % Minimum number of cores to request for calculation
JobSettings.Mempernode = '0'; % Memory request for server (default = '-1', max per core = '0', eg '3G' for cedar or 3gb for sockeye)

%% Auxiliary Job Settings (select only one)
Skip_Minimization = false; % When true, always skip the minimization step
Project_Directory_Name = 'Molten_Salts_MD';
Submit_Jobs = true; % Set to true to submit MD jobs to batch script, otherwise just produce input files.
Find_Min_Params = true; % When true, finds lowest energy parameters for IC based on Data_Types. When false, uses input IC
Find_Similar_Params = true; % If no exact minimized geometry is found for the input model, find geometry for a similar model
Data_Types = 1; % Allowed data types for automatic search of initial conditions (0 = normal, 1 = cell optimized, 2 = full optimized, 3 = atom optimized only)
N_atoms = 10000; % Minimum number of atoms to include in super cell
Table_Length = 3; % How far should tabulated potentials extend in nm
Table_StepSize = 0.0005; % Step size of tabulated potentials in nm
CoordType = 'gro'; % Either pdb, gro, or g96 (use g96 for extra precision)
Trajectory_Step = 200; %[ps] Only write frames into portable trajectory (.xtc) when t MOD dt = first time
Delete_MDPout = false; % Automatically delete MDP out file if true
Delete_MDlog = false; % Delete MD log file if true
Delete_ConfOut = false; % Delete the output configuration if true
Delete_TRR = false; % Delete the trr file after used if true
Delete_TPR = false; % Delete the tpr file after used if true
Delete_Supercell = false; % Automatically delete the supercell file if true
Delete_Backups = true; % Automatically delete any gromacs backup files found if true

%% Structure modifications
% Expand a lattice parameter of the supercell by this factor.
% Creates an empty volume upon expansion. Values less than 1 not defined.
Expand_a_SC = 1.0;
Expand_b_SC = 1.0;
Expand_c_SC = 1.0;

%% Thermostat Options
Thermostat = 'nose-hoover'; % Options: 'no' 'berendsen' 'nose-hoover' 'andersen' 'andersen-massive' 'v-rescale' (set NO for NVE)
Target_T = 298.15; % Target temperature in kelvin. Does not apply when thermostat option 'no' is chosen
Time_Constant_T = 0.4; %[ps] time constant for coupling T. Should be 20*Nsttcouple*timestep
Nsttcouple = 10; %[ps] The frequency for coupling the temperature. 

%% Barostat Options
Barostat = 'Berendsen'; % Options: 'no' 'Berendsen' 'Parrinello-Rahman' 'MTTK' (set NO for NVT)
Target_P = 1.0; % Target pressure in bar
Time_Constant_P = 1; %[ps] time constant for coupling P
Nstpcouple = 10; %[ps] The frequency for coupling the pressure. 
Compressibility = 4.5e-5; % [bar^-1] The compressibility

%% Simulated Annealing Options
Annealing = 'single'; % Options: 'no' 'single' 'periodic'
Annealing_Times = [0   500  10000 20000 30000]; % [ps] A list with the number of annealing reference/control points used
Annealing_Temps = [300 1000 1000  300   300]; % [K] A list of temperatures at the annealing reference/control points used. Must be equal in length to previous line.

%% MDP OUTPUT OPTIONS
Output_Coords = 2000; % Number of steps between outputting coordinates
Output_Coords_Compressed = 0; % Number of steps between outputting coordinates in compressed format
Output_Velocity = 2000; % Number of steps between outputting velocities
Output_Forces = 0; % Number of steps between outputting forces
Calc_Energies = 100; % Number of steps that elapse between calculating the energies.
Output_Energies = 1000; % Number of steps that else between writing energies to energy file.
Update_NeighbourList = 10; % Frequency to update the neighbor list. When this is 0, the neighbor list is made only once.

%% TOPOLOGY SETTINGS
Top_gen_pairs = 'no'; % Automatically generate pairs
Top_fudgeLJ = 1.0; % Rescale LJ interaction by this amount for 1-4 bonded atoms
Top_fudgeQQ = 1.0; % Rescale Coulomb interaction by this amount for 1-4 bonded atoms

%% MDP SETTINGS
MDP_Trajectory_Time = 30; % Trajectory time in nanoseconds. Set to 0 for single point energy calculation.
MDP_dt = 0.002; % Time step in ps for md type calculations
MDP_integrator = 'md'; % What type of calculation is run for single point energy calculations (steep = energy min, md = molecular dynamics)
MDP_LJtol = 1e-5; % When doing PME for VdW-interactions, this is used to control the relative strength of the dispersion potential at rvdw in the same way as ewald-rtol controls the electrostatic potential.
MDP_CutOffScheme = 'Verlet'; % Either 'group' or 'Verlet' (does NOT apply to tabulated potentials, these are set to group)
MDP_VerletBT = 0.005; %  (0.005) [kJ mol-1 ps-1]This sets the maximum allowed error for pair interactions per particle caused by the Verlet buffer, which indirectly sets rlist unless set to -1, in which case rlist will be used.
MDP_CoulombType = 'PME'; % Define the type of coulomb potential used. One of 'PME' or 'Cut-off'
MDP_VDWType = 'Cut-off'; % Define the type of van der waals potential used. One of 'PME' or 'Cut-off'
MDP_RList_Cutoff = 1.5; % nm. This should be larger or equal to RCoulomb/RVDW
MDP_RCoulomb_Cutoff = 1.2; % nm. if set to less than 0, then Rc = a;
MDP_RVDW_Cutoff = 1.2; % nm. note that rlist ? rCoulomb = RVDW when using Verlet and VerletBT = -1
MDP_Fourier_Spacing = 0.12; % used 0.1 in minimization. Default 0.12 nm. Grid dimensions in PME are controlled with fourierspacing
MDP_PME_Order = 4; % Interpolation order for PME. 4 equals cubic interpolation (default).
MDP_Ewald_rtol = 1e-5; %1e-7 Default (1e-5) The relative strength of the Ewald-shifted direct potential at rcoulomb. Decreasing this will give a more accurate direct sum, but then you need more wave vectors for the reciprocal sum.
MDP_nsteps = MDP_Trajectory_Time*1000/MDP_dt; % Number of steps to perform before ending (should be 0 for single point energy calculations)
MDP_Initial_T = 298.15; % Initial termpature at which to generate velocities
MDP_continuation = 'no'; % Yes = do not apply constraints to the start configuration and do not reset shells, useful for exact coninuation and reruns
MDP_Num_Groups = 2;
MDP_vdw_modifier = 'Potential-shift-Verlet'; % Potential-shift-Verlet, Potential-shift, None, Force-switch, Potential-switch

%% Pre-Minimization Settings
Maintain_Symmetry = true; % maintains cell symmetry when true
OptPos = true; % Optimize both position and lattice parameters when true
OptPos_SingleStage = false; % When true, optimize positions and lattice parameters simultaneously by numerical gradients
emtol = 0.088439; %1e-3 [kJ mol-1 nm-1] The position minimization is converged (per cycle) when the maximum force is smaller than this value
MaxCycles = 1000; % Maximum number of optimization cycles.
Energy_Tol = 1; % kJ/mol. Convergence criteria 1: must have change in energy between cycles is less than this value
Gradient_Tol_RMS = 2.5; %kJ/(mol A) Convergence criteria 2: must have RMS cell parameter gradient less than this value for convergence
Gradient_Tol_Max = 5; %kJ/(mol A) Convergence criteria 3: must have maximum cell parameter gradient less than this value for convergence
Gamma_Init = 1; % Initial step size coefficient (multiplied by numerical derivative along each dimension) used for moving along gradient
Gamma_Multiplier = 2; %1.5; % Scale Gamma by this value after each cycle (should be greater than 1)
Max_Step_size = 0.5; % Maximum possible step size in angstroms
alpha = 1/2;
beta_JC = 0.1; % Reduce Gamma by this multiple when step size is too large (for JC models).
beta_TF = 0.1; % Reduce Gamma by this multiple when step size is too large (for TF models).
MaxLineTries = 10; % Maximum number of tries to compute lower energy point before decreasing numerical derivative step size
Max_cycle_restarts = 10; % Maximum number of cycle restarts
E_Unphys = -Inf; % unphysical energy cutoff
DelE_Unphys = Inf; % Unphysical gradient cutoff
MinMDP.nsteps_point = 0; % Number of steps to perform before ending (should be 0 for single point energy calculations)
MinMDP.point_integrator = 'md'; % What type of calculation is run for single point energy calculations (steep = energy min, md = molecular dynamics)
MinMDP.dt = 0.002; % Time step in ps for md type calculations
MinMDP.min_integrator = 'steep'; % 'steep', 'cg', or 'l-bfgs'
MinMDP.nsteps_min = 1000; % Number of steps to perform before stopping for energy minimization runs
MinMDP.CutOffScheme = 'Verlet'; % Either 'group' or 'Verlet' (does NOT apply to tabulated potentials, these are set to group)
MinMDP.VerletBT = -1; % This sets the maximum allowed error for pair interactions per particle caused by the Verlet buffer, which indirectly sets rlist unless set to -1, in which case rlist will be used.

%% Begin Code
Ang_per_nm = 10; % Angstroms per nm
Longest_Cutoff = max([MDP_RList_Cutoff MDP_RCoulomb_Cutoff MDP_RVDW_Cutoff]);

% Generate bash script for batch job
[Batch_Template,qsub_cmd,gmx,Server] = MD_Batch_Template(JobSettings);

if ispc % for testing
    passlog = ' ^&^> ';
    Maindir = ['C:\Users\Hayden\Documents\Patey_Lab\' Project_Directory_Name];
    home = ['C:\Users\Hayden\Documents\Patey_Lab\' LiX_Directory]; % windows laptop
    gmx = 'wsl source ~/.bashrc; gmx_d';
    sys = @(inp) system(inp); 
elseif isunix
    passlog = ' &> ';
    if strcmpi(Server,'ced') || strcmpi(Server,'cdr') || strcmpi(Server,'sea') || strcmpi(Server,'pod')
        Maindir = ['/home/' CC_Username '/project/' Project_Directory_Name];
        home = ['/home/' CC_Username filesep LiX_Directory]; % Cedar/Graham/orcinus
        sys = @(inp) system(inp);
        cd(['/home/' CC_Username '/project']);
    elseif ~isempty(regexp(Server,'se[0-9]','ONCE')) || strcmpi(Server,'log')
        Maindir = ['/home/' CWL_Username '/scratch/' Project_Directory_Name];
        home = ['/home/' CWL_Username filesep LiX_Directory]; % Sockeye
        sys = @(inp) system(inp);
        cd(['/home/' CWL_Username '/scratch']);
    elseif strcmpi(Server,'bel')
        Maindir = ['/home/' CC_Username '/project/' Project_Directory_Name];
        home = ['/home/' CC_Username filesep LiX_Directory]; % Beluga
        sys = @(inp) system_def(inp); % Needed to circumvent error
        cd(['/home/' CC_Username '/project']);
    elseif strcmpi(Server,'pat')
        Maindir = ['/media/user/project/' Project_Directory_Name];
        home = ['/home/user/' LiX_Directory]; % Lab PC
        sys = @(inp) system(inp); 
    end
else
    error('Unknown machine type.')
end

if Potential_Input
    % Input Error checks
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
    
    if ~ischar(ModelID) && ~isnumeric(ModelID)
        error('Input variable ''ModelID'' must be a character or numeric value')
    elseif isnumeric(ModelID)
        ModelID = num2str(ModelID);
    end
    
    N = 1;
else
    % Crystal types
    if ~iscell(Structures)
        Structures = {Structures};
    end

    % Salt/ion types
    if ~iscell(Salts)
        Salts = {Salts};
    end

    % Determine number of Salts, Structures, Models, Dispersion Scaling,
    % Epsilon Scaling, Thermostat, Barostat
    N_Salt = length(Salts);
    N_Struc = length(Structures);
    N_Model = length(Models);
    N_Damp = length(Damping_Funcs);
    N_Param = length(TF_Paramset);
    N_Disp = length(Scale_Dispersion);
    N_Rep = length(Scale_Repulsion);
    N_MMD = length(Scale_MM_Dispersion);
    N_MXD = length(Scale_MX_Dispersion);
    N_XXD = length(Scale_XX_Dispersion);
    N_Epsil = length(Scale_Epsilon);
    N_Sigma = length(Scale_Sigma);
    N_Alpha = length(Scale_Alpha);
    
    N = N_Salt*N_Struc*N_Model*N_Damp*N_Param*N_Disp*N_Rep*...
        N_MMD*N_MXD*N_XXD*N_Epsil*N_Sigma*N_Alpha; % Total number of jobs

    IDX = combvec(1:N_Salt,1:N_Struc,1:N_Model,1:N_Damp,1:N_Param,1:N_Disp,1:N_Rep,...
        1:N_MMD,1:N_MXD,1:N_XXD,1:N_Epsil,1:N_Sigma,1:N_Alpha); % Vector of all possible job combinations
end

% Load topology template location
Topology_Template = fullfile(home,'templates','Gromacs_Templates',...
'Topology.template');

% Load Topology template
Topology_Template = fileread(Topology_Template);

% Add in global parameters to Topology template
Topology_Template = strrep(Topology_Template,'##GENPAIRS##',Top_gen_pairs);
Topology_Template = strrep(Topology_Template,'##FUDGELJ##',num2str(Top_fudgeLJ));
Topology_Template = strrep(Topology_Template,'##FUDGEQQ##',num2str(Top_fudgeQQ));

% Load mdp template location
MDP_Template = fullfile(home,'templates','Gromacs_Templates',...
'MDP_MD.template');

% Load MDP template
MDP_Template = fileread(MDP_Template);

% Add in global parameters to MDP template
MDP_Template = strrep(MDP_Template,'##NSTEPS##',pad(num2str(MDP_nsteps),18));
MDP_Template = strrep(MDP_Template,'##INTEGR##',pad(MDP_integrator,18));
MDP_Template = strrep(MDP_Template,'##TIMEST##',pad(num2str(MDP_dt),18));
MDP_Template = strrep(MDP_Template,'##CONTINUE##',pad(MDP_continuation,18));
MDP_Template = strrep(MDP_Template,'##REFTINIT##',pad(num2str(MDP_Initial_T),18));
MDP_Template = strrep(MDP_Template,'##COULOMB##',pad(MDP_CoulombType,18));
MDP_Template = strrep(MDP_Template,'##FOURIER##',pad(num2str(MDP_Fourier_Spacing),18));
MDP_Template = strrep(MDP_Template,'##PMEORDER##',pad(num2str(MDP_PME_Order),18));
MDP_Template = strrep(MDP_Template,'##EWALDTOL##',pad(num2str(MDP_Ewald_rtol),18));
MDP_Template = strrep(MDP_Template,'##LISTUPDATE##',pad(num2str(Update_NeighbourList),18));
MDP_Template = strrep(MDP_Template,'##POSOUT##',pad(num2str(Output_Coords),18));
MDP_Template = strrep(MDP_Template,'##POSOUTCOMP##',pad(num2str(Output_Coords_Compressed),18));
MDP_Template = strrep(MDP_Template,'##VELOUT##',pad(num2str(Output_Velocity),18));
MDP_Template = strrep(MDP_Template,'##FORCEOUT##',pad(num2str(Output_Forces),18));
MDP_Template = strrep(MDP_Template,'##ENOUT##',pad(num2str(Output_Energies),18));
MDP_Template = strrep(MDP_Template,'##CALCE##',pad(num2str(Calc_Energies),18));

if ~Submit_Jobs
    disp('Jobs will not sumit to batch system.')
	disp('To submit batch jobs automatically, set Submit_Job vaiable to true.')
end

% Ensemble type
if strcmpi(Thermostat,'no')
    if strcmpi(Barostat,'no')
        Ensemble = 'NVE';
    else
        Ensemble = 'NPH';
    end
else
    if strcmpi(Barostat,'no')
        Ensemble = 'NVT';
    else
        Ensemble = 'NPT';
    end
end

% Current Date
Date = char(datetime('today','Format','yyyy-MM-dd'));

% Loop through Job List
for idx = 1:N
    
    %% Generating Basic Job Parameters
    if Potential_Input
        Table_Req = any(Parameters(:) <= 0) || Col_Par ~= 2 || ...
            boolean(CRDamping) || boolean(C6Damping);
        
        % Model name
        Model_Scaled = ModelName(Model,Damp,0,1,1,1,1,1,1,1,1);
        Model_Scaled = [Model_Scaled '_' ModelID];
    else
        Salt = Salts{IDX(1,idx)};
        Structure = Structures{IDX(2,idx)};
        Model = Models{IDX(3,idx)};
        Damp = Damping_Funcs(IDX(4,idx));
        TF_Param = TF_Paramset(IDX(5,idx));
        S_D = Scale_Dispersion(IDX(6,idx));
        S_R = Scale_Repulsion(IDX(7,idx));
        S_MMD = Scale_MM_Dispersion(IDX(8,idx));
        S_MXD = Scale_MX_Dispersion(IDX(9,idx));
        S_XXD = Scale_XX_Dispersion(IDX(10,idx));
        S_E = Scale_Epsilon(IDX(11,idx));
        S_S = Scale_Sigma(IDX(12,idx));
        S_A = Scale_Alpha(IDX(13,idx));
        Scaling_Params = [S_D S_R S_E S_S S_MMD S_XXD S_MXD S_A];

        % Boolean: check if parameters can be input into JC without a table;
        Table_Req = (S_D <= 0) || (S_R <= 0) || ~ismembertol(1.0,S_MMD,1e-5) ...
            || ~ismembertol(1.0,S_XXD,1e-5) || ~ismembertol(1.0,S_MXD,1e-5) ...
            || Damp ~= 0;
        
        % Generate name for model with current scaling/damping parameters
        Model_Scaled = ModelName(Model,Damp,TF_Param,S_D,S_R,S_MMD,S_MXD,...
            S_XXD,S_E,S_S,S_A);
    end
    
    % Load Default Geometry info for structure
    Geometry = Default_IC(Structure);
    
    % Get Metal and Halide info from Current Salt
    [Metal,Halide] = Separate_Metal_Halide(Salt);
    Metal_Info = elements('Sym',Metal);
    Halide_Info = elements('Sym',Halide);
    
    if isempty(Model_Scaled)
        continue
    end
    
    % Find minimum lattice parameter for this salt/structure/model (or use initial ones)
    if Find_Min_Params
        [Geometry,Found_DataMatch] = FindGeometry(Geometry,Salt,Structure,...
            Model_Scaled,home,Data_Types,Find_Similar_Params);
    else
        Found_DataMatch = false;
    end
    if Skip_Minimization
        Found_DataMatch = true;
    end
    
    % Update Directory
    JobName = [Date '_' Geometry.Label '_' Model_Scaled '_' Ensemble];
    TaskName = [Geometry.Label '_' Model_Scaled '_' Ensemble];
    WorkDir = fullfile(Maindir,Salt,JobName);
    
    % Create directory if it does not exist
    if ~exist(WorkDir,'dir')
        mkdir(WorkDir)
    end
    
    % Copy Topology Template
    Topology_Text = Topology_Template;
    
    % Copy MDP Template
    MDP_Text = MDP_Template;
    
    % Copy batch job Template
    Batch_Text = Batch_Template;
    
    % Structure Template filename
    Coordinate_File = fullfile(home,'templates',[upper(CoordType) '_Templates'],...
        [Structure '.' CoordType]);

    % Load Template text
    Coordinate_Text = fileread(Coordinate_File);

    %% Adding Ions into templates
    % Insert element info into topology template
    Topology_Text = strrep(Topology_Text,'##MET##',pad(Metal,2));
    Topology_Text = strrep(Topology_Text,'##METZ##',pad(num2str(Metal_Info.atomic_number),3));
    Topology_Text = strrep(Topology_Text,'##METMASS##',pad(num2str(Metal_Info.atomic_mass),7));
    
    Topology_Text = strrep(Topology_Text,'##HAL##',pad(Halide,2));
    Topology_Text = strrep(Topology_Text,'##HALZ##',pad(num2str(Halide_Info.atomic_number),3));
    Topology_Text = strrep(Topology_Text,'##HALMASS##',pad(num2str(Halide_Info.atomic_mass),7));
    
    % Insert salt components into MDP template
    MDP_Text = strrep(MDP_Text,'##MET##',pad(Metal,3));
    MDP_Text = strrep(MDP_Text,'##HAL##',pad(Halide,3));
        
    % Add Metal and Halide symbols
    if strcmp(CoordType,'gro')
        Met = pad(Metal,2,'left');
        Hal = pad(Halide,2,'left');
    elseif strcmp(CoordType,'g96')
        Met = pad(Metal,2,'right');
        Hal = pad(Halide,2,'right');
    end
    Coordinate_Text = strrep(Coordinate_Text,'##MET##',Met);
    Coordinate_Text = strrep(Coordinate_Text,'##HAL##',Hal);
        
    % Calculate size of supercell
    N_Supercell = ceil((N_atoms/Geometry.N)^(1/3));
    
    if strcmp(Model,'TF')

        % Define the function type as 1 (required for custom functions)
        Topology_Text = strrep(Topology_Text,'##NBFUNC##','1');

        % Define the combination rules (Lorenz-berthelot)
        Topology_Text = strrep(Topology_Text,'##COMBR##','1');

        % Define all the parameters as 1.0 (already included in potentials)
        Topology_Text = strrep(Topology_Text,'##METMETC##',pad('1.0',10));
        Topology_Text = strrep(Topology_Text,'##HALHALC##',pad('1.0',10));
        Topology_Text = strrep(Topology_Text,'##METHALC##',pad('1.0',10));
        Topology_Text = strrep(Topology_Text,'##METMETA##','1.0');
        Topology_Text = strrep(Topology_Text,'##HALHALA##','1.0');
        Topology_Text = strrep(Topology_Text,'##METHALA##','1.0');

        % Generate tables of the TF potential
        if Potential_Input
            [TF_U_MX, TF_U_MM, TF_U_XX] = TF_Potential_Generator(0,...
                Table_Length,Table_StepSize,Salt,Parameters,false,...
                MDP_vdw_modifier,MDP_RVDW_Cutoff,CRDamping,Damp);
        else
            [TF_U_MX, TF_U_MM, TF_U_XX] = TFd_Potential_Generator(0,Table_Length,...
                Table_StepSize,Salt,false,Scaling_Params,MDP_vdw_modifier,MDP_RVDW_Cutoff,...
                Damp,TF_Paramset,Geometry,WorkDir);
        end

        TableName = [JobName '_Table'];
        TableFile_MX = fullfile(WorkDir,[TableName '.xvg']);
        TableFile_MM = fullfile(WorkDir,[TableName '_' Metal '_' Metal '.xvg']);
        TableFile_XX = fullfile(WorkDir,[TableName '_' Halide '_' Halide '.xvg']);

        % Save tables into current directory
        fidMX = fopen(TableFile_MX,'wt');
        fwrite(fidMX,regexprep(TF_U_MX,'\r',''));
        fclose(fidMX);

        fidMM = fopen(TableFile_MM,'wt');
        fwrite(fidMM,regexprep(TF_U_MM,'\r',''));
        fclose(fidMM);

        fidXX = fopen(TableFile_XX,'wt');
        fwrite(fidXX,regexprep(TF_U_XX,'\r',''));
        fclose(fidXX);

        % Modify the MDP file
        MDP_Text = strrep(MDP_Text,'##VDWTYPE##',pad('user',18));
        MDP_Text = strrep(MDP_Text,'##CUTOFF##',pad('group',18));
        MDP_Text = regexprep(MDP_Text,'ewald-rtol-lj.+?\n','');
        MDP_Text = regexprep(MDP_Text,'lj-pme-comb-rule.+?\n','');
        MDP_Text = regexprep(MDP_Text,'verlet-buffer-tolerance.+?\n','');
        MDP_RList_Setup = MDP_RCoulomb_Cutoff;
        MDP_RCoulomb_Setup = MDP_RCoulomb_Cutoff;
        MDP_RVDW_Setup = MDP_RVDW_Cutoff;
        MDP_Text = strrep(MDP_Text,'##RLIST##',pad(num2str(MDP_RList_Setup),18));
        MDP_Text = strrep(MDP_Text,'##RCOULOMB##',pad(num2str(MDP_RCoulomb_Setup),18));
        MDP_Text = strrep(MDP_Text,'##RVDW##',pad(num2str(MDP_RVDW_Setup),18));

        % OpenMP does not work with group cutoff
        Batch_Text = strrep(Batch_Text,'$OMP_NUM_THREADS','1');
        mpiproc = regexp(Batch_Text,'mpiprocs=([0-9]+)','tokens','ONCE');
        ompthreads = regexp(Batch_Text,'ompthreads=([0-9]+)','tokens','ONCE');
        if ~isempty(mpiproc) && ~isempty(ompthreads)
            CR = str2double(mpiproc{1})*str2double(ompthreads{1});
            Batch_Text = strrep(Batch_Text,['mpiprocs=' mpiproc{1}],['mpiprocs=' num2str(CR)]);
            Batch_Text = strrep(Batch_Text,['ompthreads=' ompthreads{1}],'ompthreads=1');
        end

    elseif contains(Model,'JC') && ~Table_Req
        switch Model
            case 'JC'
                WaterModel = 'SPC/E';
            case 'JC3P'
                WaterModel = 'TIP3P';
            case 'JC4P'
                WaterModel = 'TIP4PEW';
        end

        TableFile_MX = '';

        % Definte the function type as 1 (LJ)
        Topology_Text = strrep(Topology_Text,'##NBFUNC##','1');

        % Define the combination rules (Lorenz-berthelot in sigma-epsilon form)
        Topology_Text = strrep(Topology_Text,'##COMBR##','2');

        % Get JC parameters
        if Potential_Input
            Met_JC_Param.sigma = Parameters(1,1);
            Hal_JC_Param.sigma = Parameters(1,2);
            Met_JC_Param.epsilon = Parameters(2,1);
            Hal_JC_Param.epsilon = Parameters(2,2);

        else
            [Met_JC_Param,Hal_JC_Param] = JC_Potential_Parameters(Metal,Halide,...
                WaterModel,false,Scaling_Params);
        end
        
        % Cross terms using Lorenz-Berthelot combining rules
        Sigma_ij = (1/2)*(Met_JC_Param.sigma + Hal_JC_Param.sigma);
        epsilon_ij = sqrt(Met_JC_Param.epsilon*Hal_JC_Param.epsilon);

        % Add parameters to topology text
        Topology_Text = strrep(Topology_Text,'##METMETC##',pad(num2str(Met_JC_Param.sigma,'%10.8f'),10));
        Topology_Text = strrep(Topology_Text,'##HALHALC##',pad(num2str(Hal_JC_Param.sigma,'%10.8f'),10));
        Topology_Text = strrep(Topology_Text,'##METHALC##',pad(num2str(Sigma_ij,'%10.8f'),10));
        Topology_Text = strrep(Topology_Text,'##METMETA##',num2str(Met_JC_Param.epsilon,'%10.8f'));
        Topology_Text = strrep(Topology_Text,'##HALHALA##',num2str(Hal_JC_Param.epsilon,'%10.8f'));
        Topology_Text = strrep(Topology_Text,'##METHALA##',num2str(epsilon_ij,'%10.8f'));

        % Modify the MDP file
        MDP_Text = strrep(MDP_Text,'##VDWTYPE##',pad(MDP_VDWType,18));
        MDP_Text = strrep(MDP_Text,'##CUTOFF##',pad(MDP_CutOffScheme,18));
        MDP_Text = regexprep(MDP_Text,'energygrp-table.+?\n','');
        MDP_Text = regexprep(MDP_Text,'ewald-rtol-lj.+?\n','');
        MDP_Text = regexprep(MDP_Text,'lj-pme-comb-rule.+?\n','');
        MDP_RList_Setup = MDP_RList_Cutoff;
        MDP_RCoulomb_Setup = MDP_RCoulomb_Cutoff;
        MDP_RVDW_Setup = MDP_RVDW_Cutoff;
        MDP_Text = strrep(MDP_Text,'##RLIST##',pad(num2str(MDP_RList_Setup),18));
        MDP_Text = strrep(MDP_Text,'##RCOULOMB##',pad(num2str(MDP_RCoulomb_Setup),18));
        MDP_Text = strrep(MDP_Text,'##RVDW##',pad(num2str(MDP_RVDW_Setup),18));

        % Add in Verlet Settings
        if strcmp(MDP_CutOffScheme,'Verlet')
            MDP_Text = strrep(MDP_Text,'##VerletBT##',pad(num2str(MDP_VerletBT),18));
        else
            MDP_Text = regexprep(MDP_Text,'verlet-buffer-tolerance.+?\n','');
        end

        % Energy conversion setting
        EnergySetting = '1 2 3 4 28 29 30 31 32 33 0';

    elseif contains(Model,'JC') && Table_Req
        switch Model
            case 'JC'
                WaterModel = 'SPC/E';
            case 'JC3P'
                WaterModel = 'TIP3P';
            case 'JC4P'
                WaterModel = 'TIP4PEW';
        end

        % Define the function type as 1 (needed for custom functions)
        Topology_Text = strrep(Topology_Text,'##NBFUNC##','1');

        % Define the combination rules (Lorenz-berthelot)
        Topology_Text = strrep(Topology_Text,'##COMBR##','1');

        % Define all the parameters as 1.0 (already included in potentials)
        Topology_Text = strrep(Topology_Text,'##METMETC##',pad('1.0',10));
        Topology_Text = strrep(Topology_Text,'##HALHALC##',pad('1.0',10));
        Topology_Text = strrep(Topology_Text,'##METHALC##',pad('1.0',10));
        Topology_Text = strrep(Topology_Text,'##METMETA##','1.0');
        Topology_Text = strrep(Topology_Text,'##HALHALA##','1.0');
        Topology_Text = strrep(Topology_Text,'##METHALA##','1.0');

        % Generate tables of the JC potential
        if Potential_Input
            [JC_U_MX, JC_U_MM, JC_U_XX] = JC_Potential_Generator(0,Table_Length,...
                Table_StepSize,Salt,Parameters,false,...
                MDP_vdw_modifier,MDP_RVDW_Cutoff,CRDamping,Damp);  
        else
            [JC_U_MX, JC_U_MM, JC_U_XX] = JCd_Potential_Generator(0,Table_Length,...
                Table_StepSize,Salt,WaterModel,false,Scaling_Params,MDP_vdw_modifier,...
                MDP_RVDW_Cutoff,Damp);
        end

        TableName = [JobName '_Table'];
        TableFile_MX = fullfile(WorkDir,[TableName '.xvg']);
        TableFile_MM = fullfile(WorkDir,[TableName '_' Metal '_' Metal '.xvg']);
        TableFile_XX = fullfile(WorkDir,[TableName '_' Halide '_' Halide '.xvg']);
        
        % Save tables into current directory
        fidMX = fopen(TableFile_MX,'wt');
        fwrite(fidMX,regexprep(JC_U_MX,'\r',''));
        fclose(fidMX);

        fidMM = fopen(TableFile_MM,'wt');
        fwrite(fidMM,regexprep(JC_U_MM,'\r',''));
        fclose(fidMM);

        fidXX = fopen(TableFile_XX,'wt');
        fwrite(fidXX,regexprep(JC_U_XX,'\r',''));
        fclose(fidXX);

        % Modify the MDP file
        MDP_Text = strrep(MDP_Text,'##VDWTYPE##',pad('user',18));
        MDP_Text = strrep(MDP_Text,'##CUTOFF##',pad('group',18));
        MDP_Text = regexprep(MDP_Text,'ewald-rtol-lj.+?\n','');
        MDP_Text = regexprep(MDP_Text,'lj-pme-comb-rule.+?\n','');
        MDP_Text = regexprep(MDP_Text,'verlet-buffer-tolerance.+?\n','');
        MDP_RList_Setup = MDP_RCoulomb_Cutoff;
        MDP_RCoulomb_Setup = MDP_RCoulomb_Cutoff;
        MDP_RVDW_Setup = MDP_RVDW_Cutoff;
        MDP_Text = strrep(MDP_Text,'##RLIST##',pad(num2str(MDP_RList_Setup),18));
        MDP_Text = strrep(MDP_Text,'##RCOULOMB##',pad(num2str(MDP_RCoulomb_Setup),18));
        MDP_Text = strrep(MDP_Text,'##RVDW##',pad(num2str(MDP_RVDW_Setup),18));

        % OpenMP does not work with group cutoff
        Batch_Text = strrep(Batch_Text,'$OMP_NUM_THREADS','1');
        mpiproc = regexp(Batch_Text,'mpiprocs=([0-9]+)','tokens','ONCE');
        ompthreads = regexp(Batch_Text,'ompthreads=([0-9]+)','tokens','ONCE');
        if ~isempty(mpiproc) && ~isempty(ompthreads)
            CR = str2double(mpiproc{1})*str2double(ompthreads{1});
            Batch_Text = strrep(Batch_Text,['mpiprocs=' mpiproc{1}],['mpiprocs=' num2str(CR)]);
            Batch_Text = strrep(Batch_Text,['ompthreads=' ompthreads{1}],'ompthreads=1');
        end

        % Energy conversion setting
        EnergySetting = '1 2 3 4 28 29 30 31 32 33 0';
    else
        disp(['Warning: Unknown model type: "' Model '", skipping model...'])
        rmdir(WorkDir)
        continue
    end

    % Update MDP file with cutoff stuff
    MDP_Text = strrep(MDP_Text,'##VDWMOD##',pad(MDP_vdw_modifier,18));
    
    % Prep thermostat stuff
    if strcmpi(Thermostat,'no')
        % Remove thermostat inputs from MDP
        MDP_Text = strrep(MDP_Text,'##THERMOSTAT##',pad(Thermostat,18));
        MDP_Text = regexprep(MDP_Text,'tc-grps +.+?\n','');
        MDP_Text = regexprep(MDP_Text,'tau-t +.+?\n','');
        MDP_Text = regexprep(MDP_Text,'nsttcouple +.+?\n','');
        MDP_Text = regexprep(MDP_Text,'ref-t +.+?\n','');

    else
        % Prepare MDP inputs
        TC = '';
        TT = '';
        for x=1:MDP_Num_Groups
            TC = [num2str(Time_Constant_T) ' ' TC];
            TT = [num2str(Target_T) ' ' TT];
        end
        TC = pad(TC,18);
        TT = pad(TT,18);

        % Update MDP template
        MDP_Text = strrep(MDP_Text,'##THERMOSTAT##',pad(Thermostat,18));
        MDP_Text = strrep(MDP_Text,'##TTIMECONST##',TC);
        MDP_Text = strrep(MDP_Text,'##NSTTCOUPLE##',pad(num2str(Nsttcouple),18));
        MDP_Text = strrep(MDP_Text,'##REFT##',TT);
    end

    % Prep Barostat stuff
    if strcmpi(Barostat,'no')
        % Remove barostat inputs from MDP
        MDP_Text = strrep(MDP_Text,'##BAROSTAT##',pad(Barostat,18));
        MDP_Text = regexprep(MDP_Text,'pcoupltype +.+?\n','');
        MDP_Text = regexprep(MDP_Text,'tau-p +.+?\n','');
        MDP_Text = regexprep(MDP_Text,'compressibility +.+?\n','');
        MDP_Text = regexprep(MDP_Text,'ref-p +.+?\n','');
        MDP_Text = regexprep(MDP_Text,'nstpcouple +.+?\n','');

    else
        % Update MDP template
        MDP_Text = strrep(MDP_Text,'##BAROSTAT##',pad(Barostat,18));
        MDP_Text = strrep(MDP_Text,'##PTIMECONST##',pad(num2str(Time_Constant_P),18));
        MDP_Text = strrep(MDP_Text,'##COMPRESS##',pad(num2str(Compressibility),18));
        MDP_Text = strrep(MDP_Text,'##NSTPCOUPLE##',pad(num2str(Nstpcouple),18));
        MDP_Text = strrep(MDP_Text,'##REFP##',pad(num2str(Target_P),18));
    end

    if strcmpi(Annealing,'no')

        % Remove annealing inputs from MDP
        MDP_Text = strrep(MDP_Text,'##ANNEALING##',pad(Annealing,18));
        MDP_Text = regexprep(MDP_Text,'annealing-npoints +.+?\n','');
        MDP_Text = regexprep(MDP_Text,'annealing-time +.+?\n','');
        MDP_Text = regexprep(MDP_Text,'annealing-temp +.+?\n','');

        % Generate MDP title
        MDP_Title = [Structure ' ' Salt ' with ' Model_Scaled ' model, ' Thermostat ...
            ' thermostat, and ' Barostat ' barostat.'];

    else                   
        if length(Annealing_Times) > 1
            AT = num2str(Annealing_Temps(1));
            At = num2str(Annealing_Times(1));
            for x = 2:length(Annealing_Times)
                AT = [AT '_' num2str(Annealing_Temps(x))];
                At = [At '_' num2str(Annealing_Times(x))];
            end
        else
            AT = num2str(Annealing_Temps);
            At = num2str(Annealing_Times);
        end

        % Prepare MDP inputs
        Anneal_points = pad(repmat([num2str(length(Annealing_Times)) ' '],1,MDP_Num_Groups),18);
        Anneal_Temps = pad(repmat([strrep(AT,'_',' ') ' '],1,MDP_Num_Groups),18);
        Anneal_Times = pad(repmat([strrep(At,'_',' ') ' '],1,MDP_Num_Groups),19);

        % Update MDP template
        MDP_Text = strrep(MDP_Text,'##ANNEALING##',pad([Annealing ' ' Annealing],18));
        MDP_Text = strrep(MDP_Text,'##ANNEALPNTS##',Anneal_points);
        MDP_Text = strrep(MDP_Text,'##ANNEALTIMES##',Anneal_Times);
        MDP_Text = strrep(MDP_Text,'##ANNEALTEMPS##',Anneal_Temps);

        % Generate MDP title
        MDP_Title = [Structure ' ' Salt ' with ' Model_Scaled ' model, ' Thermostat ...
            ' thermostat, and ' Barostat ' barostat.'];
    end

    % Replace title in MDP file
    MDP_Text = strrep(MDP_Text,'##TITLE##',MDP_Title);
    
    %% Prepare remaining input files for job

    % Topology filename and directory
    Topology_File = fullfile(WorkDir,[JobName '.top']);

    % Save MDP file in current directory
    MDP_File = fullfile(WorkDir,[JobName '.mdp']);
    fidMDP = fopen(MDP_File,'wt');
    fwrite(fidMDP,regexprep(MDP_Text,'\r',''));
    fclose(fidMDP);
    
    % Unit Cell Filename
    UnitCellFile = fullfile(WorkDir,[JobName '_UnitCell.' CoordType]);

    % Supercell Filename
    SuperCellFile = fullfile(WorkDir,[JobName '.' CoordType]);
    
    % Convert to N x N x N supercell
    N_Cell = Geometry.N;
    N = num2str(N_Supercell);
    N_total = (N_Supercell^3)*N_Cell;
    
    % Create name for mdpout file
    MDPout_File = fullfile(WorkDir,[JobName '_out.mdp']);

    % Grompp log file
    GrompLog_File = fullfile(WorkDir,[JobName '_Grompplog.log']);

    % Prepare trajectory file
    Trajectory_File = fullfile(WorkDir,[JobName '.tpr']);
    
    % Add coordinates in xyz space (lattice parameter-dependent)
    if Found_DataMatch
        
        % Save number of atoms into .mat file
        NumberFile = fullfile(WorkDir,[JobName '.mat']);
        save(NumberFile,'N_total','N_Cell')
        
        % Add number of unit cells to topology file
        Topology_Text = strrep(Topology_Text,'##N##',num2str(N_Supercell));
        
        Coordinate_Text = AddCartesianCoord(Coordinate_Text,Geometry,1,false,CoordType);
        Batch_Text = strrep(Batch_Text,['##PREMIN##' newline],'');
        
        % Save unit cell .gro file into current directory
        fid = fopen(UnitCellFile,'wt');
        fwrite(fid,regexprep(Coordinate_Text,'\r',''));
        fclose(fid);

        Supercell_command = [gmx ' genconf -f ' windows2unix(UnitCellFile) ...
             ' -o ' windows2unix(SuperCellFile) ' -nbox ' N ' ' N ' ' N];
        [errcode,output] = sys(Supercell_command);

        if errcode ~= 0
            disp(output);
            error(['Error creating supercell with genconf. Problem command: ' newline Supercell_command]);
        end
        
        % Expand supercell by requested amount
        a_sc = num2str(Expand_a_SC*Geometry.a*N_Supercell/10,'%10.8f'); % supercell a length in nm
        b_sc = num2str(Expand_b_SC*Geometry.b*N_Supercell/10,'%10.8f'); % supercell b length in nm
        c_sc = num2str(Expand_c_SC*Geometry.c*N_Supercell/10,'%10.8f'); % supercell c length in nm
        
        bc = num2str(Geometry.alpha,'%10.4f');
        ac = num2str(Geometry.beta,'%10.4f');
        ab = num2str(Geometry.gamma,'%10.4f');
        
        Expand_command = [gmx ' editconf -f ' windows2unix(SuperCellFile) ...
             ' -o ' windows2unix(SuperCellFile) ' -box ' a_sc ' ' b_sc ' ' c_sc ' ' ...
             '-noc -angles ' bc ' ' ac ' ' ab];
        [errcode,output] = sys(Expand_command);
        
        if errcode ~= 0
            disp(output);
            error(['Error expanding supercell with editconf. Problem command: ' newline Expand_command]);
        end
        
        % Generate topology file
        Atomlist = copy_atom_order(SuperCellFile);
        Topology_Text = strrep(Topology_Text,'##LATOMS##',Atomlist);
        fidTOP = fopen(Topology_File,'wt');
        fwrite(fidTOP,regexprep(Topology_Text,'\r',''));
        fclose(fidTOP);
        
        GROMPP_command = [gmx ' grompp -c ' windows2unix(SuperCellFile) ...
            ' -f ' windows2unix(MDP_File) ' -p ' windows2unix(Topology_File) ...
            ' -o ' windows2unix(Trajectory_File) ' -po ' windows2unix(MDPout_File) ...
            ' -maxwarn 1' passlog windows2unix(GrompLog_File)];
        [errcode,~] = sys(GROMPP_command);
        
        % Catch error in grompp
        if errcode ~= 0
            error(['Error running GROMPP. Problem command: ' newline GROMPP_command]);
        end
        
    else
        % Save number of atoms into .mat file
        MinDir = fullfile(WorkDir,'Minimization');
        if ~exist(MinDir,'dir')
            mkdir(MinDir);
        end
        JobFile = fullfile(MinDir,'TempJobInfo.mat');
        save(JobFile);
        Batch_Text = strrep(Batch_Text,'##PREMIN##',['matlab -batch "MD_Preminimization(''' MinDir ''')"']);
    end

    % Prepare mdrun command
    Log_File = fullfile(WorkDir,[JobName '.log']);

    Energy_file = fullfile(WorkDir,[JobName '.edr']);

    TRR_File = fullfile(WorkDir,[JobName '.trr']);

    ConfOut_File = fullfile(WorkDir,[JobName 'OutConf.' CoordType]);

    CheckPoint_File = fullfile(WorkDir,[JobName '.cpt']);

    XTC_File = fullfile(WorkDir,[JobName '.xtc']);

    mdrun_command = [gmx ' mdrun -s ' windows2unix(Trajectory_File) ...
        ' -o ' windows2unix(TRR_File) ' -g ' windows2unix(Log_File) ...
        ' -e ' windows2unix(Energy_file) ' -c ' windows2unix(ConfOut_File) ...
        ' -cpo ' windows2unix(CheckPoint_File)];

    if contains(Model_Scaled,'TF') || Table_Req % TF potential requires table
        mdrun_command = [mdrun_command ' -table ' windows2unix(TableFile_MX)];
    end

    % Remove files if requested
    % Clean up MDP out file if set  
    cleanup_command = '';
    if Delete_MDPout
        cleanup_command = ['rm ' MDPout_File newline]; %#ok<*UNRCH>
    end
    if Delete_Supercell
        cleanup_command = [cleanup_command 'rm ' SuperCellFile newline];
    end
    if Delete_MDlog
        cleanup_command = [cleanup_command 'rm ' Log_File newline];
    end
    if Delete_ConfOut
        cleanup_command = [cleanup_command 'rm ' ConfOut_File newline];
    end
    if Delete_TRR
        cleanup_command = [cleanup_command 'rm ' TRR_File newline];
    end
    if Delete_TPR
        cleanup_command = [cleanup_command 'rm ' Trajectory_File newline];
    end
    if Delete_Backups
        cleanup_command = [cleanup_command 'find ' WorkDir ...
            ' -name ''#*#'' -delete' newline]; %#ok<*AGROW>
    end
    
    % After the run, convert into compressed trajectory file
    cleanup_command = [cleanup_command ...
        gmx ' trjconv -f ' TRR_File ' -s ' Trajectory_File ...
        ' -pbc atom -ur tric -o ' XTC_File ' -dt ' num2str(Trajectory_Step)];

    % Place into batch script
    Batch_Text = strrep(Batch_Text,'##MDRUN##',mdrun_command);
    Batch_Text = strrep(Batch_Text,'##CLEANUP##',cleanup_command);
    Batch_Text = strrep(Batch_Text,'##TASKNAME##',TaskName);
    Batch_Text = strrep(Batch_Text,'##ERROR##',[WorkDir filesep TaskName]);
    Batch_Text = strrep(Batch_Text,'##DIRECTORY##',WorkDir);

    % Open and save batch script
    fidBS = fopen(fullfile(WorkDir,[JobName '.subm']),'wt');
    fwrite(fidBS,regexprep(Batch_Text,'\r',''));
    fclose(fidBS);

    % Submit job
    disp('Job input files produced for:')
    disp([Salt ' ' TaskName])
    if ~ispc && Submit_Jobs
        sys([qsub_cmd ' ' fullfile(WorkDir,[JobName '.subm'])]);
        disp('Job submitted.')
    end
end

end