% MATLAB Parallel wrapper for Structure_Minimization to run minimization 
% of all structures in parallel
function Output_Array =  Structure_Minimization_Par(Salt,Model,Parameters,OptPos)
    %Structures = {'Rocksalt' 'Wurtzite' 'Sphalerite' 'CsCl' 'NiAs' 'BetaBeO' 'FiveFive'};
    Structures = {'Rocksalt' 'CsCl'};
    N = length(Structures);
    Output_Array = cell(1,N);

    Parcores = feature('numcores');
    if isempty(gcp('nocreate'))
        parpool(min(Parcores,N));
    end

    parfor idx = 1:N
        Structure = Structures{idx};
        filename = [Salt '_' Structure '_' Model '.log'];
        diary(filename);
        Output_Array{idx} = Structure_Minimization(Salt,Structure,Model,Parameters,OptPos);
    end
end

