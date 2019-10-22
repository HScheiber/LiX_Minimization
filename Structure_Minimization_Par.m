% MATLAB Parallel wrapper for Structure_Minimization to run minimization 
% of all structures in parallel
function Output_Array =  Structure_Minimization_Par(Salt,Model,Parameters,OptPos)
    % Pre-construct output array and list of structures to compute
    Structures = {'Rocksalt' 'Wurtzite' 'Sphalerite' 'CsCl' 'NiAs' 'BetaBeO' 'FiveFive'};
    
    N = length(Structures);
    Data_Array = cell(1,N);
    
    % Create a unique filename based on the current time for the log file
    filename = [char(datetime(now,'ConvertFrom','datenum','format','yyyy-MM-dd_HH.mm.ss')) ...
        '_' Salt '_' Model '.log'];

    % Set up matlab parallel features
    Parcores = feature('numcores');
    if ~isempty(gcp('nocreate'))
        delete(gcp);
    end
    ppool = parpool(min(Parcores,N));
    partitions = 1:N;
    f(partitions) = parallel.FevalFuture;
    
    % Run in parallel with parfavel
    for idx = partitions
        f(idx) = parfeval(ppool,@Structure_Minimization,1,Salt,Structures{idx},Model,Parameters,OptPos);
    end
    wait(f);
    
    % Collect outputs and save diaries to single string variable
    Diary_Text = '';
    for idx = partitions
        Data_Array{idx} = f(idx).OutputArguments;
        Diary_Text = [Diary_Text f(idx).Diary newline newline]; %#ok<AGROW>
    end
    
    % Save diaries to text file    
    fid = fopen(filename,'wt');
    fprintf(fid,Diary_Text');
    fclose(fid);
    
    % Determine if wurtzite has converged to 5-5
    Wurtz_ind = strcmp(Structures,'Wurtzite');
    if sum(Wurtz_ind) > 0
        Wurtzite_Data = Data_Array{Wurtz_ind}{1};
        delta_z = mod(Wurtzite_Data(7) - Wurtzite_Data(10),1);
        if abs(delta_z-0.5) < 1e-2
            Structures{Wurtz_ind} = 'FiveFive';
        end
    end
    
    % Determmine if BetaBeO has moved to rocksalt
    bBeO_ind = strcmp(Structures,'BetaBeO');
    if sum(bBeO_ind) > 0
        BetaBeO_Data = Data_Array{bBeO_ind}{1};
        a_minus_c = abs(BetaBeO_Data(2) - BetaBeO_Data(4)) < 1e-2;
        x1 = mod(BetaBeO_Data(5),1) - (1/4) < 1e-2;
        x2 = mod(BetaBeO_Data(8),1) - (1/4) < 1e-2;
        y1 = mod(BetaBeO_Data(6),1) - (3/4) < 1e-2;
        y2 = mod(BetaBeO_Data(9),1) - (3/4) < 1e-2;
        z1 = mod(BetaBeO_Data(7),1) < 1e-2;
        z2 = mod(BetaBeO_Data(10),1) < 1e-2;
        if a_minus_c && x1 && x2 && y1 && y2 && z1 && z2
            Structures{bBeO_ind} = 'Rocksalt';
        end
    end
    
    % Find lowest energy structure
    Energies = nan(1,N);
    for i = 1:N
        Energies(i) = Data_Array{i}{1}(1);
    end
    [~,Ind] = min(Energies);
    MinStruct = Structures{Ind};
    
    % Output data only if lowest energy structure is rocksalt
    if strcmp(MinStruct,'Rocksalt')
        Output_Array = Data_Array{Ind}{1};
    else
        Output_Array = nan(1,10);
    end
end