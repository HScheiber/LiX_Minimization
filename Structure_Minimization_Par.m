% MATLAB Parallel wrapper for Structure_Minimization to run minimization 
% of all structures in parallel
function Output_Array =  Structure_Minimization_Par(Salt,Model,Parameters,OptPos)
    Structures = {'Rocksalt' 'Wurtzite' 'Sphalerite' 'CsCl' 'NiAs' 'BetaBeO' 'FiveFive'};
    N = length(Structures);
    Output_Array = cell(1,N);

    Parcores = feature('numcores');
    if isempty(gcp('nocreate'))
        parpool(Parcores);
    end

    parfor idx = 1:N
        Structure = Structures{idx};
        Output_Array{idx} = Structure_Minimization(Salt,Structure,Model,Parameters,OptPos);
    end
end

