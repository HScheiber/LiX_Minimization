% MATLAB Parallel wrapper for Structure_Minimization to run minimization 
% of all structures in parallel
function Output_Array =  Structure_Minimization_Par(Salt,Model,Parameters,OptPos)
    % Pre-construct output array and list of structures to compute
    Structures = {'Rocksalt' 'Wurtzite' 'Sphalerite' 'CsCl' 'NiAs' 'BetaBeO' 'FiveFive'};
    N = length(Structures);
    Output_Array = cell(1,N);
    
    % Create a unique filename based on the current time for the log file
    filename = [char(datetime(now,'ConvertFrom','datenum','format','yyyy-MM-dd-HH:mm:ss')) ...
        '_' Salt '_' Model '_All.log'];

    % Set up matlab parallel features
    Parcores = feature('numcores');
    if isempty(gcp('nocreate'))
        ppool = parpool(min(Parcores,N));
    else
        ppool = gcp;
    end
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
        Output_Array{idx} = f(idx).OutputArguments;
        Diary_Text = [Diary_Text f(idx).Diary newline newline];
    end
    
    % Save diaries to text file    
    fid = fopen(filename,'wt');
    fprintf(fid,Diary_Text');
    fclose(fid);
    
end