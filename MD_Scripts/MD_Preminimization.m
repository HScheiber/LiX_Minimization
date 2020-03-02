%% INFO ABOUT INPUTS
% Models is a string variable. One of 'JC' or 'TF'. Sets the mathematical form of the model to use.
%
% Salts is either a string or cell array. Can be any of: 'LiF' 'LiCl' 'LiBr' 'LiI' 'NaCl'
%
% Structures is either a string or cell array. Can be any of: 'Rocksalt' 'Wurtzite' 'Sphalerite' 'CsCl' 'NiAs' 'BetaBeO' 'FiveFive'
%
% OptPos is a Boolean. If true, then the optimization algorithm optimizes
% for lattice parameters AND fractional coordinates. If false, the
% algorithm only optimizes for lattice parameters.
%
% CRDamping is a Boolean, short for Close-Range Damping. If true, this adds
% a damping function for the potential at close range. NOTE that using this
% will force the potential to be defined by a table in GROMACS, which runs
% the software somewhat slower than otherwise.

% C6Damping is an integer. Adds medium range dispersion-only damping:
% BJ damping: https://www.chemie.uni-bonn.de/pctc/mulliken-center/software/dft-d3/man.pdf
% Other damping functions: https://www.pamoc.it/kw_dfd.html
% 0 = no (default) damping
% 1 = BJ/rational damping (same as in D3(BJ), damps to a constant. Fairly
% weak damping)
% 2 = Tang Damping (Mid strength damping, damps to zero at r=0)
% 3 = MMDRE Damping function (very weak damping)
% 4 = PAMoC Damping function (weak damping)
% 5 = EHFSK Damping function (strong damping)
% 6 = WY damping function (strongest damping)
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
function MD_Preminimization(Directory)

%% Move to input directory and load input variables
cd(Directory)
IP = load('TempJobInfo.mat');

%% Save Minimization Diary
DiaryFile = fullfile(IP.WorkDir,[IP.JobName '.minlog']);
diary(DiaryFile);

%% Begin Code
Longest_Cutoff = max([IP.MDP_RList_Cutoff IP.MDP_RCoulomb_Cutoff IP.MDP_RVDW_Cutoff]);

% Load mdp template location
MDP_Template = fullfile(IP.home,'templates','Gromacs_Templates',...
'MDP.template');

% Load MDP template
MDP_Template = fileread(MDP_Template);

% Add in global parameters to MDP template
MDP_Template = strrep(MDP_Template,'##NSTEPS##',pad(num2str(IP.MinMDP.nsteps_point),18));
MDP_Template = strrep(MDP_Template,'##INTEGR##',pad(IP.MinMDP.point_integrator,18));
MDP_Template = strrep(MDP_Template,'##TIMEST##',pad(num2str(IP.MinMDP.dt),18));

% Determine type of optimization
if IP.OptPos
    IP.OptTxt = 'FULLOPT';
else
    IP.OptTxt = 'CELLOPT';
end

% Insert salt components into MDP template
MDP_Template = strrep(MDP_Template,'##MET##',IP.Metal);
MDP_Template = strrep(MDP_Template,'##HAL##',IP.Halide);

% Select symmetry settings
if IP.Maintain_Symmetry
    switch IP.Structure
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
    IP.OptTxt = [IP.OptTxt '_SG1'];
    DOF = {'a' 'b' 'c' 'alpha' 'beta' 'gamma'};
    DOF_txt = {'a' 'b' 'c' char(945) char(946) char(947)};
    DOF_Units = {[' ' char(0197)] [' ' char(0197)] [' ' char(0197)] ...
        [' ' char(0176)] [' ' char(0176)] [' ' char(0176)]};
end

% Extra DOF for Optimization of Fractional Coordinates in 1 step
FC_DOFs = {[IP.Metal '_x'] [IP.Metal '_y'] [IP.Metal '_z'] [IP.Halide '_x'] [IP.Halide '_y'] [IP.Halide '_z']};
if IP.OptPos_SingleStage && IP.OptPos
    DOF(end+1:end+6) = FC_DOFs;
    DOF_txt(end+1:end+6) = FC_DOFs;
    DOF_Units(end+1:end+6) = {'' '' '' '' '' ''};
end

% How many Degrees of freedom in gradient
N_DOF = length(DOF);

% Get Structure Template filename
Coordinate_File = fullfile(IP.home,'templates',...
    [upper(IP.CoordType) '_Templates'],...
    [IP.Structure '.' IP.CoordType]);

% Load Template text from disk
Coordinate_text = fileread(Coordinate_File);

% Add Metal and Halide symbols
if strcmp(IP.CoordType,'gro')
    Met = pad(IP.Metal,2,'left');
    Hal = pad(IP.Halide,2,'left');
elseif strcmp(IP.CoordType,'g96')
    Met = pad(IP.Metal,2,'right');
    Hal = pad(IP.Halide,2,'right');
end
Coordinate_text = strrep(Coordinate_text,'##MET##',Met);
Coordinate_text = strrep(Coordinate_text,'##HAL##',Hal);

% Update File Base Name
FileBase = [IP.Salt '_' IP.Geometry.Label '_' IP.Model_Scaled '_' IP.OptTxt];

% Update Topology and MDP files
MDP_Template = strrep(MDP_Template,'##COULOMB##',pad('PME',18));
MDP_Template = strrep(MDP_Template,'##FOURIER##',pad(num2str(IP.MDP_Fourier_Spacing),18));
MDP_Template = strrep(MDP_Template,'##PMEORDER##',pad(num2str(IP.MDP_PME_Order),18));
MDP_Template = strrep(MDP_Template,'##EWALDTOL##',pad(num2str(IP.MDP_Ewald_rtol),18));

if strcmpi(IP.Model,'TF')

    beta = IP.beta_TF;
    
    % Modify the MDP file
    MDP_Template = strrep(MDP_Template,'##VDWTYPE##',pad('user',18));
    MDP_Template = strrep(MDP_Template,'##VDWMOD##',pad(IP.MDP_vdw_modifier,18));
    MDP_Template = strrep(MDP_Template,'##CUTOFF##',pad('group',18));
    MDP_Template = regexprep(MDP_Template,'ewald-rtol-lj.+?\n','');
    MDP_Template = regexprep(MDP_Template,'lj-pme-comb-rule.+?\n','');
    MDP_Template = regexprep(MDP_Template,'verlet-buffer-tolerance.+?\n','');
    
    % Energy conversion setting
    EnergySetting = '1 2 3 4 28 29 30 31 32 33 0';
    
elseif contains(IP.Model,'JC') && ~IP.Table_Req
    beta = IP.beta_JC;

    % Modify the MDP file
    MDP_Template = strrep(MDP_Template,'##VDWTYPE##',pad(IP.MDP_VDWType,18));
    MDP_Template = strrep(MDP_Template,'##VDWMOD##',pad(IP.MDP_vdw_modifier,18));
    MDP_Template = strrep(MDP_Template,'##CUTOFF##',pad(IP.MDP_CutOffScheme,18));
    MDP_Template = regexprep(MDP_Template,'energygrp-table.+?\n','');
    MDP_Template = regexprep(MDP_Template,'ewald-rtol-lj.+?\n','');
    MDP_Template = regexprep(MDP_Template,'lj-pme-comb-rule.+?\n','');

    % Add in Verlet Settings
    if strcmp(IP.MDP_CutOffScheme,'Verlet')
        MDP_Template = strrep(MDP_Template,'##VerletBT##',pad(num2str(IP.MDP_VerletBT),18));
    else
        MDP_Template = regexprep(MDP_Template,'verlet-buffer-tolerance.+?\n','');
    end

    % Energy conversion setting
    EnergySetting = '1 2 3 4 28 29 30 31 32 33 0';

elseif contains(IP.Model,'JC') && IP.Table_Req

    beta = IP.beta_JC;

    % Modify the MDP file
    MDP_Template = strrep(MDP_Template,'##VDWTYPE##',pad('user',18));
    MDP_Template = strrep(MDP_Template,'##VDWMOD##',pad(IP.MDP_vdw_modifier,18));
    MDP_Template = strrep(MDP_Template,'##CUTOFF##',pad('group',18));
    MDP_Template = regexprep(MDP_Template,'ewald-rtol-lj.+?\n','');
    MDP_Template = regexprep(MDP_Template,'lj-pme-comb-rule.+?\n','');
    MDP_Template = regexprep(MDP_Template,'verlet-buffer-tolerance.+?\n','');

    % Energy conversion setting
    EnergySetting = '1 2 3 4 28 29 30 31 32 33 0';
end


if ispc
    DirectoryUnix = windows2unix(Directory);
    rm_command = ['wsl find ' DirectoryUnix ' -iname "#*#" -delete'];
    del_command = ['wsl rm -r ' DirectoryUnix '/*'];
else
    rm_command = ['find ' Directory ' -iname "#*#" -delete'];
    del_command = ['rm -r ' Directory '/*'];
end

TotalTimer = tic;
disp(['Beginning ' IP.Salt ' ' IP.Structure ' ' IP.Model_Scaled ' Optimization...'])
disp('********************** Convergence Requirements ***********************')
disp(['Max ' char(916) 'E between cycles: ' num2str(IP.Energy_Tol,'%4.4E') ' a.u.']);
disp(['Max RMS ' char(8711) 'E: ' num2str(IP.Gradient_Tol_RMS,'%4.4E') ' a.u. / ' char(0197)]);
disp(['Max component ' char(8711) 'E: ' num2str(IP.Gradient_Tol_Max,'%4.4E') ' a.u. / ' char(0197)]);

%% Begin Optimization Loop
% Loop broken when convergence criteria is met or max cycles reached
skip_results = false;
auto_h = true;
Cycle_restarts = 0;
Gamma = IP.Gamma_Init;
for Index = 1:IP.MaxCycles
    Restart_cycle = false;
    % Calculate reasonable step size
    if auto_h
        h = zeros(1,N_DOF);

        IP.N_Supercell = MDOptimizationLoop(IP,IP.Geometry,Directory,MDP_Template,...
            Longest_Cutoff,Coordinate_text,EnergySetting);
        E = GrabEnergy(Directory,FileBase);

        disp(['********************Cycle ' num2str(Index) ' Initial Conditions********************']);
        for Didx = 1:N_DOF
            if ~ismember(DOF{Didx},FC_DOFs)
                h(Didx) = nuderst(IP.Geometry.(DOF{Didx}));
                disp(['Lattice Parameter ' DOF_txt{Didx} ' = ' ...
                    num2str(IP.Geometry.(DOF{Didx}),'%2.8f') DOF_Units{Didx} '.' ...
                    ' Num. Der. Step Size: ' char(948) '(' DOF_txt{Didx} ') = ' num2str(h(Didx),'%2.8E') DOF_Units{Didx} '.']);
            else
                h(Didx) = nuderst(1);
            end
        end
        disp(['Fractional Coordinates for ' IP.Metal ': ']);
        disp(num2str(IP.Geometry.FC_Metal(:,:),'%2.8f '))
        disp(['Fractional Coordinates for ' IP.Halide ': ']);
        disp(num2str(IP.Geometry.FC_Halide(:,:),'%2.8f '))
        disp(['Initial E = ' num2str(E,'%4.10f') ' kJ/mol']);
        disp(['Step size coefficient: ' num2str(Gamma,'%4.4E') ])
        Telap = datestr(seconds(toc(TotalTimer)),'HH:MM:SS');
        disp(['Time Elapsed: ' Telap])
        disp('******************************************************************')
    else
        auto_h = true;
    end

    % Check for unphysical energy
    if E < IP.E_Unphys
        disp(['Warning: Unphysical energy detected after ' num2str(Index) ' cycles.'])
        disp('Model may have local minimum at complete overlap of opposite charges.')
        disp('Otherwise poor initial conditions may cause this.')
        disp('Removing output files.')
        system(del_command);
        skip_results = true;
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
        GeoMinus = IP.Geometry;

        if ismember(DOF{Didx},FC_DOFs)
            % For FC degrees of freedom
            MX = DOF{Didx}(1:2);
            xyz = strcmp(DOF{Didx}(end),{'x' 'y' 'z'});
            if strcmp(MX,pad(IP.Metal,2,'_'))
                GeoMinus.FC_Metal(1,xyz) = mod(GeoMinus.FC_Metal(1,xyz) - delta,1);
            else
                GeoMinus.FC_Halide(1,xyz) = mod(GeoMinus.FC_Halide(1,xyz) - delta,1);
            end

            % Grab asymmetric unit FC
            AsymFC_Metal = GeoMinus.FC_Metal(1,:);
            AsymFC_Halide = GeoMinus.FC_Halide(1,:);

            % Reseat fractional coordinates outside of asymmetric unit
            [GeoMinus.FC_Metal,GeoMinus.FC_Halide] = ...
                UnitCell_FractionalCoords(AsymFC_Metal,AsymFC_Halide,Structure);
        else
            % For lattice parameter degrees of freedom
            GeoMinus.(CurDOF) = GeoMinus.(CurDOF) - delta;

            % Update transformation matrix
            GeoMinus.Transform = GenTransformMatrix(GeoMinus);

            % Maintain symmetry if set
            if IP.Maintain_Symmetry && strcmp(CurDOF,'a')
                [GeoMinus.b,GeoMinus.c] = SymmetryAdapt(GeoMinus.a,...
                    GeoMinus.b,GeoMinus.c,IP.Structure); %#ok<*UNRCH>
            end

        end

        % Generate energy at new point
        MDOptimizationLoop(IP,GeoMinus,Directory,MDP_Template,...
            Longest_Cutoff,Coordinate_text,EnergySetting);
        E_Minus = GrabEnergy(Directory,FileBase);

        Telap = datestr(seconds(toc(TotalTimer)),'HH:MM:SS');
        disp(['Finite step ' CurDOFtxt ' -' num2str(delta,'%2.2E') ...
            CurDOFunit ' Calculated. ' char(916) 'E(' CurDOFtxt '-' char(948) ...
            ') = ' num2str(E_Minus-E,'%+4.4E') ' kJ/mol. Time Elapsed: ' Telap]);

        % Get energy at plus step
        GeoPlus = IP.Geometry;

        if ismember(DOF{Didx},FC_DOFs)
            % For FC degrees of freedom
            MX = DOF{Didx}(1:2);
            xyz = strcmp(DOF{Didx}(end),{'x' 'y' 'z'});
            if strcmp(MX,pad(IP.Metal,2,'_'))
                GeoPlus.FC_Metal(1,xyz) = mod(GeoPlus.FC_Metal(1,xyz) + delta,1);
            else
                GeoPlus.FC_Halide(1,xyz) = mod(GeoPlus.FC_Halide(1,xyz) + delta,1);
            end

            % Grab asymmetric unit FC
            AsymFC_Metal = GeoPlus.FC_Metal(1,:);
            AsymFC_Halide = GeoPlus.FC_Halide(1,:);

            % Reseat fractional coordinates outside of asymmetric unit
            [GeoPlus.FC_Metal,GeoPlus.FC_Halide] = ...
                UnitCell_FractionalCoords(AsymFC_Metal,AsymFC_Halide,IP.Structure);
        else
            % For lattice parameter degrees of freedom
            GeoPlus.(CurDOF) = GeoPlus.(CurDOF) + delta;

            % Update transformation matrix
            GeoPlus.Transform = GenTransformMatrix(GeoPlus);

            % Maintain symmetry if set
            if IP.Maintain_Symmetry && strcmp(CurDOF,'a')
                [GeoPlus.b,GeoPlus.c] = SymmetryAdapt(GeoPlus.a,...
                    GeoPlus.b,GeoPlus.c,IP.Structure); %#ok<*UNRCH>
            end
        end

        % Generate Energy at new point
        MDOptimizationLoop(IP,GeoPlus,Directory,MDP_Template,...
            Longest_Cutoff,Coordinate_text,EnergySetting);
        E_Plus = GrabEnergy(Directory,FileBase);

        Telap = datestr(seconds(toc(TotalTimer)),'HH:MM:SS');
        disp(['Finite step ' CurDOFtxt ' +' num2str(delta,'%2.2E') ...
            CurDOFunit ' Calculated. ' char(916) 'E(' CurDOFtxt '+' char(948) ...
            ') = ' num2str(E_Plus-E,'%+4.4E') ' kJ/mol. Time Elapsed: ' Telap]);

        % 3-Point Numerical derivative wrt current DOF
        Gradient(Didx) = (1/(2*delta))*(-E_Minus + E_Plus);
        system(rm_command);
    end
    disp([char(8711) 'E = ' num2str(Gradient,'%5.4E  ') ' kJ/mol ' char(0197)])

    if max(Gradient) > IP.DelE_Unphys
        disp(['Warning: Unphysical energy gradient detected after ' num2str(Index) ' cycles.'])
        disp('Model may have local minimum at complete overlap of opposite charges.')
        disp('Otherwise poor initial conditions may cause this.')
        disp('Removing output files.')
        system(del_command);
        skip_results = true;
        break
        % Check for convergence on cycle 1
    elseif Index == 1 && (rms(Gradient) < IP.Gradient_Tol_RMS) && ...
            (max(abs(Gradient)) < IP.Gradient_Tol_Max) && ~IP.OptPos
        % If gradient convergence criteria are met, end loop
        disp('Energy convergence reached on cycle 1 gradient test.')
        system(rm_command);
        break
    end

    % Move one step in direction of steepest descent
    for StepInd = 1:IP.MaxLineTries

        GeomNew = IP.Geometry;

        for Didx = 1:N_DOF

            CurDOF = DOF{Didx};
            CurGrad = Gradient(Didx);

            % Step size for this degree of freedom along gradient
            % component
            ShiftStep = -sign(CurGrad)*min(abs(Gamma*CurGrad),IP.Max_Step_size);

            if ismember(CurDOF,FC_DOFs)
                % For FC degrees of freedom
                MX = CurDOF(1:2);
                xyz = strcmp(CurDOF(end),{'x' 'y' 'z'});
                if strcmp(MX,pad(IP.Metal,2,'_'))
                    GeomNew.FC_Metal(1,xyz) = mod(GeomNew.FC_Metal(1,xyz) + ShiftStep,1);
                else
                    GeomNew.FC_Halide(1,xyz) = mod(GeomNew.FC_Halide(1,xyz) + ShiftStep,1);
                end
            else
                % For lattice parameter degrees of freedom
                GeomNew.(CurDOF) = IP.Geometry.(CurDOF) + ShiftStep;
            end
        end

        if IP.OptPos_SingleStage
            % Grab asymmetric unit FC
            AsymFC_Metal = GeomNew.FC_Metal(1,:);
            AsymFC_Halide = GeomNew.FC_Halide(1,:);

            % Reseat fractional coordinates outside of asymmetric unit
            [GeomNew.FC_Metal,GeomNew.FC_Halide] = ...
                UnitCell_FractionalCoords(AsymFC_Metal,AsymFC_Halide,IP.Structure);
        end

        if IP.Maintain_Symmetry
            [GeomNew.b,GeomNew.c] = SymmetryAdapt(GeomNew.a,...
                GeomNew.b,GeomNew.c,IP.Structure);
        end

        % Update transformation matrix
        GeomNew.Transform = GenTransformMatrix(GeomNew);

        % Recalculate energy at new point
        IP.N_Supercell = MDOptimizationLoop(IP,GeomNew,Directory,MDP_Template,...
            Longest_Cutoff,Coordinate_text,EnergySetting);
        E_New = GrabEnergy(Directory,FileBase);
        E_New_LP = E_New;

        if E_New > E - min(IP.alpha*Gamma*norm(Gradient)^2,50)
            Gamma = beta*Gamma;

            Telap = datestr(seconds(toc(TotalTimer)),'HH:MM:SS');
            disp(['Optimizing step size... Time Elapsed: ' Telap])
            [~,~] = system(rm_command);
        else
            Telap = datestr(seconds(toc(TotalTimer)),'HH:MM:SS');
            disp(['Finite step in direction of steepest descent calculated. Time Elapsed: ' Telap]);
            disp([char(916) 'E = ' num2str(E_New-E,'%+4.4E') ' kJ/mol']);
            disp(['RMS ' char(8711) 'E = ' num2str(rms(Gradient),'%4.4E') ' kJ/mol ' char(0197)])
            disp(['|Max(' char(8711) 'E)| = ' num2str(max(abs(Gradient)),'%4.4E') ' kJ/mol ' char(0197)])
            break
        end

        if StepInd >= IP.MaxLineTries
            h = max(h./10,eps*1000);
            Gamma = IP.Gamma_Init;
            Restart_cycle = true;
            system(rm_command);
            break
        end
    end

    if Restart_cycle
        Cycle_restarts = Cycle_restarts+1;
        if Cycle_restarts >= IP.Max_cycle_restarts
            disp(['Unable to find lower energy point after ' num2str(Max_cycle_restarts) ' attempts.'])
            disp('Stopping here.')
            if (rms(Gradient) < IP.Gradient_Tol_RMS) && (max(abs(Gradient)) < IP.Gradient_Tol_Max)
                disp('Gradient convergence criteria fulfilled.')
                disp(['Energy convergence reached after ' num2str(Index) ' cycles.'])
            else
                disp('Gradient convergence criteria NOT fulfilled.')
                disp(['Search haulted after ' num2str(Index) ' cycles.'])
                disp('Check initial conditions. Poor initial conditions may cause this.')
                disp('Removing output files.')
                system(del_command);
                skip_results = true;
            end

            break
        else
            auto_h = false;
            disp('Unable to find a point of lower energy in chosen direction, starting next cycle with modified numerical derivative step size.');
            continue
        end
    end

    % Check for unphysical energy
    if (E_New < IP.E_Unphys) || (max(Gradient) > IP.DelE_Unphys)
        disp(['Warning: Unphysical energy detected after ' num2str(Index) ' cycles.'])
        disp('Model may have local minimum at complete overlap of opposite charges.')
        disp('Otherwise poor initial conditions may cause this.')
        disp('Removing output files.')
        system(del_command);
        skip_results = true;
        break
    end

    % Optimize wrt positions
    if IP.OptPos && ~IP.OptPos_SingleStage
        % Save MDP variable for later
        Temp_ModelOld = MDP_Template;

        MDP_Template = strrep(MDP_Template,'= md                ; What type of calculation is run',...
            ['= ' IP.MinMDP.min_integrator newline 'emtol                    = ' num2str(IP.emtol)]);
        MDP_Template = regexprep(MDP_Template,'nsteps                   = [0-9|\.|\-]+',...
            ['nsteps                   = ' num2str(IP.MinMDP.nsteps_min)]);

        % Rerun energy calc and break loop
        EnergySettingAlt = '1 2 3 4 25 26 27 28 29 30 0'; 
        N_Supercell_out = MDOptimizationLoopFC(IP,GeomNew,Directory,MDP_Template,...
            Longest_Cutoff,Coordinate_text,EnergySettingAlt);
        
        E_Old = E_New;
        E_New = GrabEnergyFinal(Directory,FileBase);
        disp(['Geometry Optimized W.R.T. Atomic Positions. ' char(916) 'E = ' num2str(E_New-E,'%+4.4E') ' kJ/mol']);
        Telap = datestr(seconds(toc(TotalTimer)),'HH:MM:SS');
        disp(['Time Elapsed: ' Telap])

        % Get new unit cell coordinates

        % Get output config filename
        OutConf_file = fullfile(Directory,[FileBase 'OutConf.' IP.CoordType]);

        % Open the unit cell structure file
        SuperCellText = fileread(OutConf_file);
        SuperCellText = regexprep(SuperCellText,'\r','');
        
        % g96 file format
        if strcmpi(IP.CoordType,'g96')
            
            Positions = regexp(SuperCellText,'(?<=POSITION\n)(.+?)(?=(\n)END)','match','ONCE');
            Coords_Data = textscan(Positions,'%*6c%*6c%*6c%*6c%15.9f%15.9f%15.9f\n',...
                GeomNew.N,'Delimiter','','whitespace','');
            
            % Get end cell parameter text
            boxtemp = regexp(SuperCellText,'(?<=BOX\n)(.+?)(?=(\n)END)','match');
            boxcoords = textscan(boxtemp{1},'%f %f %f %f %f %f %f %f %f',...
                'Delimiter',' ','MultipleDelimsAsOne',true);
            
            % Check for empty coords
            indx = cellfun('isempty',boxcoords); % true for empty cells
            boxcoords(indx) = {0}; % replace by a cell with a zero
        % gro file format
        elseif strcmpi(IP.CoordType,'gro')

            Coords_Data = textscan(SuperCellText,'%*5d%*-5s%*5s%*5d%8.3f%8.3f%8.3f%*8.4f%*8.4f%*8.4f\n',...
                GeomNew.N,'Delimiter','','whitespace','','HeaderLines',2);
            
            % Get end cell parameter text
            SkipLines = (N_Supercell_out^3)*GeomNew.N + 2;
            boxcoords = textscan(SuperCellText,'%f %f %f %f %f %f %f %f %f',...
                'Delimiter',' ','MultipleDelimsAsOne',true,'HeaderLines',SkipLines);
            
            % Check for empty coords
            indx = cellfun(@isnan,boxcoords); % true for empty cells
            boxcoords(indx) = {0}; % replace by a cell with a zero
        end

        % Get xyz coordinates of the unit cell
        idxX = (GeomNew.N)/2 + 1;
        Met_XYZ = [Coords_Data{1}(1:idxX-1) Coords_Data{2}(1:idxX-1) Coords_Data{3}(1:idxX-1)];
        Hal_XYZ = [Coords_Data{1}(idxX:end) Coords_Data{2}(idxX:end) Coords_Data{3}(idxX:end)];

        % Lattice vectors
        a_vec = [boxcoords{1} boxcoords{4} boxcoords{5}]/N_Supercell_out;
        b_vec = [boxcoords{6} boxcoords{2} boxcoords{7}]/N_Supercell_out;
        c_vec = [boxcoords{8} boxcoords{9} boxcoords{3}]/N_Supercell_out;

        if (strcmp(IP.Structure,'Wurtzite') || strcmp(IP.Structure,'NiAs')) && VecAngle(a_vec,b_vec) < 90.0
            b_vec = [-boxcoords{6} boxcoords{2} boxcoords{7}]/N_Supercell_out;
        end

        % Get transformation Matrix
        TM = [a_vec ; b_vec ; c_vec];
        GeomNew.Transform = [a_vec/norm(a_vec) ; b_vec/norm(b_vec) ; c_vec/norm(c_vec)];

        % Convert to Fractional coordinates
        GeomNew.FC_Metal = mod(Met_XYZ/TM,1);
        GeomNew.FC_Halide = mod(Hal_XYZ/TM,1);

        % Copy to unit cell
        Unit_cell_text = AddCartesianCoord(Coordinate_text,GeomNew,1,false,IP.CoordType);

        % Overwrite unit cell
        UnitCell_file = fullfile(Directory,[FileBase '_UnitCell.' IP.CoordType]);
        fidUC = fopen(UnitCell_file,'wt');
        fwrite(fidUC,regexprep(Unit_cell_text,'\r',''));
        fclose(fidUC);

        % Reload old MDP variable
        MDP_Template = Temp_ModelOld;
    end

    % Remove backup files
    system(rm_command);

    Telap = datestr(seconds(toc(TotalTimer)),'HH:MM:SS');
    disp('************************************')
    disp(['Cycle ' num2str(Index) ' complete.'])
    disp(['Time Elapsed: ' Telap])
    disp('************************************')

    % Check if energy converged
    if (abs(E_New - E) < IP.Energy_Tol) && ...
            (rms(Gradient) < IP.Gradient_Tol_RMS) && ...
            (max(abs(Gradient)) < IP.Gradient_Tol_Max)
        % If all convergence criteria are met, end loop
        IP.Geometry = GeomNew;
        disp(['Energy convergence reached after ' num2str(Index) ' cycles.' newline 'Final recalculation of energy...'])
        IP.N_Supercell = MDOptimizationLoop(IP,IP.Geometry,Directory,MDP_Template,...
            Longest_Cutoff,Coordinate_text,EnergySetting);
        
        E = GrabEnergy(Directory,FileBase);
        system(rm_command);
        break
    % If energy converged but not gradients
    elseif (abs(E_New_LP - E) < IP.Energy_Tol)
        Gamma = IP.Gamma_Init; % Re-initialize Gamma
        E = E_New;
        IP.Geometry = GeomNew;
    % Detect unphysical energy
    elseif (E_New < IP.E_Unphys) || (max(Gradient) > IP.DelE_Unphys)
        disp(['Warning: Unphysical energy detected after ' num2str(Index) ' cycles.'])
        disp('Model may have local minimum at complete overlap of opposite charges.')
        disp('Otherwise poor initial conditions may cause this.')
        disp('Removing output files.')
        system(del_command);
        skip_results = true;
        break
    % Otherwise a normal decrease in energy
    elseif E_New < E
        Gamma = Gamma*IP.Gamma_Multiplier;
        E = E_New;
        IP.Geometry = GeomNew;
    % If energy increases
    elseif E_New > E
        disp('Warning: Total energy increased after geometry optimization of atomic positions.')
        disp('Convergence may not have been reached by GROMACS atomic position minimization step.')
        disp('Reverting to intial fractional coordinates.');
        Gamma = Gamma*IP.Gamma_Multiplier;
        E = E_Old;
        GeomNew.FC_Metal = IP.Geometry.FC_Metal;
        GeomNew.FC_Halide = IP.Geometry.FC_Halide;
        IP.Geometry = GeomNew;
    end

    if Index < IP.MaxCycles
        disp(['Beginning Cycle ' num2str(Index+1) '.'])
    else
        disp(['Convergence NOT reached after ' num2str(Index) ' cycles. Stopping.'])
    end                        
end

Telap = datestr(seconds(toc(TotalTimer)),'HH:MM:SS');
disp(['Completed: ' IP.Salt ' ' IP.Structure ' ' IP.Model_Scaled ' Geometry Optimization. Time ' Telap])
if ~skip_results
    disp(['Final Optimized Energy is ' num2str(E,'%4.10f') ' kJ/mol'])

    for Didx = 1:N_DOF
        if ~ismember(DOF{Didx},FC_DOFs)
            disp(['Lattice Parameter ' DOF_txt{Didx} ' = ' ...
                num2str(IP.Geometry.(DOF{Didx}),'%2.8f') DOF_Units{Didx} '.']);
        end
    end
    disp(['Fractional Coordinates for ' IP.Metal ': ']);
    disp(num2str(IP.Geometry.FC_Metal(:,:),'%2.8f '))
    disp(['Fractional Coordinates for ' IP.Halide ': ']);
    disp(num2str(IP.Geometry.FC_Halide(:,:),'%2.8f '))
    disp('Components of Final Gradient:')
    disp(num2str(Gradient,'%+4.4E  '))
else
    error('Could not minimize structure.');
end
diary off

% Add number of unit cells to topology file (this will not change)
IP.Topology_Text = strrep(IP.Topology_Text,'##N##',num2str(IP.N_Supercell));

% Save number of atoms into .mat file (this wont change)
NumberFile = fullfile(IP.WorkDir,[IP.JobName '.mat']);
N_Cell = IP.N_Cell;
N_total = (IP.N_Supercell^3)*N_Cell;
save(NumberFile,'N_total','N_Cell');
N = num2str(IP.N_Supercell);

%% Geometry Editing of minimized cell plus final gromacs minimization

% Add unit cell coordinates
IP.Coordinate_Text = AddCartesianCoord(IP.Coordinate_Text,IP.Geometry,1,false,IP.CoordType);

% Save unit cell .gro file into main directory
fid = fopen(IP.UnitCellFile,'wt');
fwrite(fid,regexprep(IP.Coordinate_Text,'\r',''));
fclose(fid);

% Create supercell
Temp_SupercellFile = fullfile(Directory,[IP.TaskName '.' IP.CoordType]);
Supercell_command = [IP.gmx ' genconf -f ' windows2unix(IP.UnitCellFile) ...
     ' -o ' windows2unix(Temp_SupercellFile) ' -nbox ' N ' ' N ' ' N];
[errcode,output] = system(Supercell_command);

if errcode ~= 0
    disp(output);
    error(['Error creating supercell with genconf. Problem command: ' newline Supercell_command]);
end

% Geometry Editing: expand supercell by requested amount
a_sc = num2str(IP.Expand_a_SC*IP.Geometry.a*IP.N_Supercell/10,'%10.8f'); % supercell a length in nm
b_sc = num2str(IP.Expand_b_SC*IP.Geometry.b*IP.N_Supercell/10,'%10.8f'); % supercell b length in nm
c_sc = num2str(IP.Expand_c_SC*IP.Geometry.c*IP.N_Supercell/10,'%10.8f'); % supercell c length in nm

% Cell angles
bc = num2str(IP.Geometry.alpha,'%10.4f');
ac = num2str(IP.Geometry.beta,'%10.4f');
ab = num2str(IP.Geometry.gamma,'%10.4f');

Expand_command = [IP.gmx ' editconf -f ' windows2unix(Temp_SupercellFile) ...
     ' -o ' windows2unix(Temp_SupercellFile) ' -box ' a_sc ' ' b_sc ' ' c_sc ' ' ...
     '-noc -angles ' bc ' ' ac ' ' ab];
[errcode,output] = system(Expand_command);

if errcode ~= 0
    disp(output);
    error(['Error expanding supercell with editconf. Problem command: ' newline Expand_command]);
end

%% Gromacs final minimize after geometry editing
MDP_Template = strrep(MDP_Template,'= md                ; What type of calculation is run',...
    ['= ' IP.MinMDP.min_integrator newline 'emtol                    = ' num2str(IP.emtol)]);
MDP_Template = regexprep(MDP_Template,'nsteps                   = [0-9|\.|\-]+',...
    ['nsteps                   = ' num2str(IP.MinMDP.nsteps_min)]);

% Determine cutoff length
R_List_Cutoff = pad(num2str(IP.MDP_RList_Setup),18);
R_Coulomb_Cutoff = pad(num2str(IP.MDP_RCoulomb_Setup),18);
R_VDW_Cutoff = pad(num2str(IP.MDP_RVDW_Setup),18);
MDP_Template = strrep(MDP_Template,'##RLIST##',R_List_Cutoff);
MDP_Template = strrep(MDP_Template,'##RCOULOMB##',R_Coulomb_Cutoff);
MDP_Template = strrep(MDP_Template,'##RVDW##',R_VDW_Cutoff);

% Save MDP file in current directory
MDP_File = fullfile(Directory,[FileBase '.mdp']);
fidMDP = fopen(MDP_File,'wt');
fwrite(fidMDP,regexprep(MDP_Template,'\r',''));
fclose(fidMDP);

% Generate topology file
Topology_File = fullfile(Directory,[FileBase '.top']);
Atomlist = copy_atom_order(Temp_SupercellFile);
Topology_text_new = strrep(IP.Topology_Text,'##LATOMS##',Atomlist);
fidTOP = fopen(Topology_File,'wt');
fwrite(fidTOP,Topology_text_new);
fclose(fidTOP);

Trajectory_File = fullfile(Directory,[FileBase '.tpr']);
MDPout_File = fullfile(Directory,[FileBase '_out.mdp']);
GrompLog_File = fullfile(Directory,[FileBase '_Grompplog.log']);

FMin_Grompp = [IP.gmx ' grompp -c ' windows2unix(Temp_SupercellFile) ...
    ' -f ' windows2unix(MDP_File) ' -p ' windows2unix(Topology_File) ...
    ' -o ' windows2unix(Trajectory_File) ' -po ' windows2unix(MDPout_File) ...
    ' -maxwarn 1' IP.passlog windows2unix(GrompLog_File)];
[state,~] = system(FMin_Grompp);
% Catch error in grompp
if state ~= 0
    error(['Error running GROMPP. Problem command: ' newline FMin_Grompp]);
else
    delete(GrompLog_File)
end

% Prepare minimization mdrun command
Log_File = fullfile(Directory,[FileBase '.log']);

Energy_file = fullfile(Directory,[FileBase '.edr']);

TRR_File = fullfile(Directory,[FileBase '.trr']);

mdrun_command = [IP.gmx ' mdrun -s ' windows2unix(Trajectory_File) ...
    ' -o ' windows2unix(TRR_File) ' -g ' windows2unix(Log_File) ...
    ' -e ' windows2unix(Energy_file) ' -c ' windows2unix(IP.SuperCellFile)];
if ~isempty(IP.TableFile_MX) % TF potential requires table
    mdrun_command = [mdrun_command ' -table ' windows2unix(IP.TableFile_MX)];
end

% Final minimization
[state,mdrun_output] = system(mdrun_command);
if state ~= 0
    disp(mdrun_output);
    error(['Error running mdrun for final minimization. Problem command: ' newline mdrun_command]);
end

% Generate final topology file for molecular dynamics
Atomlist = copy_atom_order(IP.SuperCellFile);
IP.Topology_Text = strrep(IP.Topology_Text,'##LATOMS##',Atomlist);
fidTOP = fopen(IP.Topology_File,'wt');
fwrite(fidTOP,regexprep(IP.Topology_Text,'\r',''));
fclose(fidTOP);

GROMPP_command = [IP.gmx ' grompp -c ' windows2unix(IP.SuperCellFile) ...
    ' -f ' windows2unix(IP.MDP_File) ' -p ' windows2unix(IP.Topology_File) ...
    ' -o ' windows2unix(IP.Trajectory_File) ' -po ' windows2unix(IP.MDPout_File) ...
    ' -maxwarn 1' IP.passlog windows2unix(IP.GrompLog_File)];
[errcode,~] = system(GROMPP_command);

% Catch error in grompp
if errcode ~= 0
    error(['Error running GROMPP. Problem command: ' newline GROMPP_command]);
end

% Return to working directory of primary job
cd(IP.WorkDir)

% Remove minimization temporary folder
system(del_command);
rmdir(Directory);
return
end
