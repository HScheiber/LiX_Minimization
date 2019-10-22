function [a,b,c,FCmet,FChal,Updated] = FindMinLatticeParam(Cry,Salt,Structure,...
    Model,home,Data_Types,Find_Similar_Params)

% Data files
Datafilename = fullfile(home,'DATA',...
      [Salt '_' Structure '_Lattice_Energies.mat']);

% Load known empirical energies
if exist(Datafilename, 'file') == 2
    X = load(Datafilename,'Data');
    Data = X.Data;
    clearvars X
else
    % If no data available, load default and return   
    a = Cry.(Structure).a;
    b = Cry.(Structure).b;
    c = Cry.(Structure).c;
    FCmet = Cry.(Structure).FC_Metal;
    FChal = Cry.(Structure).FC_Halide;
    Updated = false;
    return
end  

%% Check if field of interest present in data structure
if ~isfield(Data,Salt)
    a = Cry.(Structure).a;
    b = Cry.(Structure).b;
    c = Cry.(Structure).c;
    FCmet = Cry.(Structure).FC_Metal;
    FChal = Cry.(Structure).FC_Halide;
    Updated = false;
    return
end
if ~isfield(Data.(Salt),Structure)
    a = Cry.(Structure).a;
    b = Cry.(Structure).b;
    c = Cry.(Structure).c;
    FCmet = Cry.(Structure).FC_Metal;
    FChal = Cry.(Structure).FC_Halide;
    Updated = false;
    return
end
if ~isfield(Data.(Salt).(Structure),Model)
    
    if Find_Similar_Params
        % Find similar model
        if contains(Model,'JC3P')
            Model = 'JC3P';
        elseif contains(Model,'JC4P')
            Model = 'JC4P';
        elseif contains(Model,'JC')
            Model = 'JC';
        elseif contains(Model,'TF')
            Model = 'TF';
        end
    end
    
     % If still no model found, use default
	if ~isfield(Data.(Salt).(Structure),Model)
        a = Cry.(Structure).a;
        b = Cry.(Structure).b;
        c = Cry.(Structure).c;
        FCmet = Cry.(Structure).FC_Metal;
        FChal = Cry.(Structure).FC_Halide;
        Updated = false;
        return
	end
end

EmpiricalData = Data.(Salt).(Structure).(Model);
a_emp = [EmpiricalData{:,1}];
b_emp = [EmpiricalData{:,2}];
c_emp = [EmpiricalData{:,3}];
E_emp = [EmpiricalData{:,7}];
FC_emp = {EmpiricalData{:,6}};
DT = [EmpiricalData{:,9}];

% Remove unselected data types
DatInd = ismember(DT,Data_Types);
a_emp = a_emp(DatInd);
b_emp = b_emp(DatInd);
c_emp = c_emp(DatInd);
E_emp = E_emp(DatInd);
FC_emp = FC_emp(DatInd);

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
    a = Cry.(Structure).a;
    b = Cry.(Structure).b;
    c = Cry.(Structure).c;
    FCmet = Cry.(Structure).FC_Metal;
    FChal = Cry.(Structure).FC_Halide;
    Updated = false;
    return
end

% Global minimum
[~,Ind] = min(E_emp);
a = a_emp(Ind);
b = b_emp(Ind);
c = c_emp(Ind);
FC = FC_emp{Ind};

% Generate fractional coordinates in primary unit cell for output
[FCmet,FChal] = UnitCell_FractionalCoords(...
    [FC{1,2} FC{1,3} FC{1,4}],...
    [FC{2,2} FC{2,3} FC{2,4}],Structure);
Updated = true;
end


