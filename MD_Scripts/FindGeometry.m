function [Geometry,Found_DataMatch] = FindGeometry(Geometry,Salt,Structure,...
    Model_Scaled,home,Data_Types,Find_Similar_Params)

Found_DataMatch = true;

% Convert structure name for model
ModelTag = strrep(Model_Scaled,'.','P');
Datafilename = fullfile(home,'data','GROMACS',...
      [Salt '_' Structure '_Lattice_Energies.mat']);

% Load known empirical energies
if exist(Datafilename, 'file') == 2
    X = load(Datafilename,'Data');
    Data = X.Data;
    clearvars X
else
    % If no data available, load default and return   
    Found_DataMatch = false;
    return
end  

%% Check if field of interest present in data structure
if ~isfield(Data,Salt)
    Found_DataMatch = false;
    return
end
if ~isfield(Data.(Salt),Structure)
    Found_DataMatch = false;
    return
end
if ~isfield(Data.(Salt).(Structure),ModelTag)
    Found_DataMatch = false;
    if Find_Similar_Params
        % Find similar model
        if contains(ModelTag,'JC3P')
            ModelTag = 'JC3P';
        elseif contains(ModelTag,'JC4P')
            ModelTag = 'JC4P';
        elseif contains(ModelTag,'JC')
            ModelTag = 'JC';
        elseif contains(ModelTag,'TF')
            ModelTag = 'TF';
        end
    end
    
     % If still no model found, use default
	if ~isfield(Data.(Salt).(Structure),ModelTag)
        return
	end
end

EmpiricalData = Data.(Salt).(Structure).(ModelTag);
a_emp = [EmpiricalData{:,1}];
b_emp = [EmpiricalData{:,2}];
c_emp = [EmpiricalData{:,3}];
E_emp = [EmpiricalData{:,7}];
FC_emp = EmpiricalData(:,6)';
DT = [EmpiricalData{:,9}];

% Remove unselected data types
DatInd = ismember(DT,Data_Types);

% No data of selected type available? Then keep all data if
% Find_Similar_Params is true
if sum(DatInd) == 0 && Find_Similar_Params
    Found_DataMatch = false;
    Data_Types = 1:4;
elseif sum(DatInd) == 0
    Found_DataMatch = false;
    return
else
    a_emp = a_emp(DatInd);
    b_emp = b_emp(DatInd);
    c_emp = c_emp(DatInd);
    E_emp = E_emp(DatInd);
    FC_emp = FC_emp(DatInd);
end

% Local minima
if ismember(0,Data_Types)
    Ind = islocalmin(E_emp);
    a_emp = a_emp(Ind);
    b_emp = b_emp(Ind);
    c_emp = c_emp(Ind);
    E_emp = E_emp(Ind);
    FC_emp = FC_emp(Ind);
end

% No local minimum found, use default settings
if isempty(a_emp)
    Found_DataMatch = false;
    return
end

% Global minimum
[~,Ind] = min(E_emp);
Geometry.a = a_emp(Ind);
Geometry.b = b_emp(Ind);
Geometry.c = c_emp(Ind);
FC = FC_emp{Ind};

% Generate fractional coordinates in primary unit cell for output
[Geometry.FC_Metal,Geometry.FC_Halide] = UnitCell_FractionalCoords(...
    [FC{1,2} FC{1,3} FC{1,4}],...
    [FC{2,2} FC{2,3} FC{2,4}],Structure);
end


