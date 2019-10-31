% MATLAB Parallel wrapper for Structure_Minimization to run minimization 
% of an input number of Parameter sets for a given Salt, Structure, and
% model in parallel.
%
%% INFO ABOUT INPUTS
% Models is a string variable. One of 'JC' or 'TF'. Sets the mathematical form of the model to use.
% Salts is either a string or cell array. Can be any of: 'LiF' 'LiCl' 'LiBr' 'LiI' 'NaCl'
% Structure is a string. Can be any of: 'Rocksalt' 'Wurtzite' 'Sphalerite' 'CsCl' 'NiAs' 'BetaBeO' 'FiveFive'
% OptPos is a Boolean. If true, then the optimization algorithm optimizes
% for lattice parameters AND fractional coordinates. If false, the
% algorithm only optimizes for lattice parameters.
%
%% Parameters is and N x M x P matrix of floats, N and M depend on the chosen model. 
% P is the number of parameter sets to run. The code will attempt to run all P
% jobs in parallel, but will only run the same number of jobs as there
% are cores on the local machine.
%
% For each 1, ..., P:
%
% M and N depend on the particular model chosen.
%
%
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
function Output_Array = Structure_Minimization_Multi(Salt,Structure,Model,Parameters,OptPos)

    % Check that Parameters is 3D matrix
    [~,~,N] = size(Parameters);
    
    % Create a unique filename based on the current time for the log file
    filename = [char(datetime(now,'ConvertFrom','datenum','format','yyyy-MM-dd_HH.mm.ss')) ...
        '_' Salt '_' Structure '_' Model '.log'];

    % Set up matlab parallel features
    Parcores = feature('numcores');
    PrefCores = min(Parcores,N);
    
    if ~isempty(gcp('nocreate'))
        Cur_Pool = gcp;
        Cur_Workers = Cur_Pool.NumWorkers;
        if Cur_Workers < PrefCores
            delete(Cur_Pool);
            ppool = parpool(min(Parcores,N));
        else
            ppool = Cur_Pool;
        end
    else
        ppool = parpool(min(Parcores,N)); 
    end
    partitions = 1:N;
    f(partitions) = parallel.FevalFuture;
    
    % Run in parallel with parfavel
    for idx = partitions
        f(idx) = parfeval(ppool,@Structure_Minimization,1,Salt,Structure,Model,Parameters(:,:,idx),OptPos);
    end
    wait(f);
    
    % Collect outputs and save diaries to single string variable
    Diary_Text = '';
    Data_Array = cell(1,N);
    for idx = partitions
        Data_Array{idx} = f(idx).OutputArguments{1};
        Diary_Text = [Diary_Text f(idx).Diary newline newline]; %#ok<AGROW>
    end
    
    % Save diaries to text file    
    fid = fopen(filename,'wt');
    fprintf(fid,Diary_Text');
    fclose(fid);
    
    % Determine if wurtzite has converged to 5-5
%     if strcmp(Structure,'Wurtzite')
%         for idx = partitions
%             Wurtzite_Data = Data_Array{idx};
%             delta_z = mod(Wurtzite_Data(7) - Wurtzite_Data(10),1);
%             if abs(delta_z-0.5) < 1e-2
%                 Data_Array{idx}(:) = nan;
%             end
%         end
%     end
    
    % Determmine if BetaBeO has moved to rocksalt
%     if strcmp(Structure,'BetaBeO')
%         for idx = partitions
%             BetaBeO_Data = Data_Array{idx};
%             a_minus_c = abs(BetaBeO_Data(2) - BetaBeO_Data(4)) < 1e-2;
%             x1 = mod(BetaBeO_Data(5),1) - (1/4) < 1e-2;
%             x2 = mod(BetaBeO_Data(8),1) - (1/4) < 1e-2;
%             y1 = mod(BetaBeO_Data(6),1) - (3/4) < 1e-2;
%             y2 = mod(BetaBeO_Data(9),1) - (3/4) < 1e-2;
%             z1 = mod(BetaBeO_Data(7),1) < 1e-2;
%             z2 = mod(BetaBeO_Data(10),1) < 1e-2;
%             if a_minus_c && x1 && x2 && y1 && y2 && z1 && z2
%                 Data_Array{idx}(:) = nan;
%             end
%         end
%     end

    Output_Array = nan(N,10);
    for idx = partitions
        Output_Array(idx,:) = Data_Array{idx};
    end

end