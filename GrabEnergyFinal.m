function Energy = GrabEnergyFinal(Directory,FileBase)

    % Find xvg file
    enfile_info = dir(fullfile(Directory,...
        [FileBase 'energies.xvg']));

    % Find total number of atoms .mat file
    totatoms_info = dir(fullfile(Directory,...
        [FileBase '.mat']));

    if isempty(enfile_info) || isempty(totatoms_info)
        error('Optimization step failed to produce required output files.');
    end

    % Get total atoms filename
    Atoms_file = fullfile(Directory,totatoms_info.name);

    % Open the mat file to get total number of atoms/molecules
    N_info = load(Atoms_file,'-mat','N_Cell','N_total');
    N_atoms = N_info.N_Cell;
    N_total_mols = N_info.N_total/2;

    % Get energy filename
    Energy_file = fullfile(Directory,enfile_info.name);

    % Import energy file as text
    Energy_text = fileread(Energy_file);

    % Find the final energy line
    Energies_Final = textscan(Energy_text, '%s', 'delimiter', '\n','MultipleDelimsAsOne',true);

    Energy_list = textscan(Energies_Final{1}{end},'%*f %f %f %f %f %f %f %f %f %f %f',...
        'Delimiter',' ','MultipleDelimsAsOne',true);
    Energy_array = cell2mat(Energy_list)./N_total_mols;

    Energy = Energy_array(4);
end