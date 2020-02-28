% Function which generates default structure parameters
function Geometry = Default_IC(Structure)

switch Structure
    case 'Rocksalt'
        Geometry.a = 5.2027;
        Geometry.b = Geometry.a;
        Geometry.c = Geometry.a;
        Geometry.FC_Metal  = [0.0 0.0 0.0];
        Geometry.FC_Halide = [1/2 1/2 1/2];
        Geometry.alpha = 90;
        Geometry.beta = 90;
        Geometry.gamma = 90;
        Geometry.Transform = eye(3);
        [Geometry.FC_Metal,Geometry.FC_Halide] = ...
            UnitCell_FractionalCoords(Geometry.FC_Metal,...
            Geometry.FC_Halide,'Rocksalt');
        Geometry.N = 8;
        Geometry.Label = 'R';
    case 'Wurtzite'
        Geometry.a = 4.0222;
        Geometry.b = Geometry.a;
        Geometry.c = sqrt(8/3)*Geometry.a; % Perfect Wurtzite c/a = sqrt(8/3);
        Geometry.FC_Metal  = [1/3 2/3 0.0]; %[1/3 2/3 0.0];
        Geometry.FC_Halide = [1/3 2/3 3/8]; %[1/3 2/3 3/8];
        Geometry.alpha = 90;
        Geometry.beta = 90;
        Geometry.gamma = 120;
        Geometry.Transform =  [1        0        0; ...
                              -sind(30) cosd(30) 0; ...
                               0        0        1];
        [Geometry.FC_Metal,Geometry.FC_Halide] = ...
            UnitCell_FractionalCoords(Geometry.FC_Metal,...
            Geometry.FC_Halide,'Wurtzite');
        Geometry.N = 4;
        Geometry.Label = 'W';
    case 'FiveFive'
        Geometry.a = 5.1815;
        Geometry.b = (3/2)*Geometry.a; % Theoretical b/a = 3/2
        Geometry.c = cosd(30)*Geometry.a; % Theoretical c/a = cosd(30)
        Geometry.FC_Metal  = [1/4 1/6 1/2];
        Geometry.FC_Halide = [1/4 1/3 0.0];
        Geometry.alpha = 90;
        Geometry.beta = 90;
        Geometry.gamma = 90;
        Geometry.Transform = eye(3);
        [Geometry.FC_Metal,Geometry.FC_Halide] = ...
            UnitCell_FractionalCoords(Geometry.FC_Metal,...
            Geometry.FC_Halide,'FiveFive');
        Geometry.N = 8;
        Geometry.Label = 'F';
    case 'NiAs'
        x = 1.391666666667; %1.7316 LiF LSDA, 1.7275 LiF PBE
        Geometry.a = 3.6318;
        Geometry.b = Geometry.a;
        Geometry.c = Geometry.a*x; 
        Geometry.FC_Metal  = [0.0 0.0 0.0];
        Geometry.FC_Halide = [1/3 2/3 1/4];
        Geometry.alpha = 90;
        Geometry.beta  = 90;
        Geometry.gamma = 120;
        Geometry.Transform =  [1        0        0; ...
                              -sind(30) cosd(30) 0; ...
                               0        0        1];
        [Geometry.FC_Metal,Geometry.FC_Halide] = ...
            UnitCell_FractionalCoords(Geometry.FC_Metal,...
            Geometry.FC_Halide,'NiAs');
        Geometry.N = 4;
        Geometry.Label = 'N';
    case 'Sphalerite'
        Geometry.a = 5.6726;
        Geometry.b = Geometry.a;
        Geometry.c = Geometry.a;
        Geometry.FC_Metal  = [0.0 0.0 0.0];
        Geometry.FC_Halide = [1/4 1/4 1/4];
        Geometry.alpha = 90;
        Geometry.beta  = 90;
        Geometry.gamma = 90;
        Geometry.Transform = eye(3);
        [Geometry.FC_Metal,Geometry.FC_Halide] = ...
            UnitCell_FractionalCoords(Geometry.FC_Metal,...
            Geometry.FC_Halide,'Sphalerite');
        Geometry.N = 8;
        Geometry.Label = 'S';
    case 'CsCl'
        Geometry.a = 3.3202;
        Geometry.b = Geometry.a;
        Geometry.c = Geometry.a;
        Geometry.FC_Metal  = [0.0 0.0 0.0];
        Geometry.FC_Halide = [1/2 1/2 1/2];
        Geometry.alpha = 90;
        Geometry.beta  = 90;
        Geometry.gamma = 90;
        Geometry.Transform = eye(3);
        [Geometry.FC_Metal,Geometry.FC_Halide] = ...
            UnitCell_FractionalCoords(Geometry.FC_Metal,...
            Geometry.FC_Halide,'CsCl');
        Geometry.N = 2;
        Geometry.Label = 'C';
    case 'BetaBeO'
        x = 0.336; %0.336 Experimental value
        y = 0.310; %0.310 Experimental value
        Geometry.a = 6.9274; % 5.8
        Geometry.b = Geometry.a;
        Geometry.c = (sqrt(3)/2)*Geometry.a; %(sqrt(3)/2)*a
        Geometry.FC_Metal  = [x 1-x 0.000];
        Geometry.FC_Halide = [y y   0.000];
        Geometry.alpha = 90;
        Geometry.beta  = 90;
        Geometry.gamma = 90;
        Geometry.Transform = eye(3);
        [Geometry.FC_Metal,Geometry.FC_Halide] = ...
            UnitCell_FractionalCoords(Geometry.FC_Metal,...
            Geometry.FC_Halide,'BetaBeO');
        Geometry.N = 8;
        Geometry.Label = 'B';
    otherwise
        error(['Unknown default structure: ' Structure])
end
Geometry.NF = Geometry.N/2;
end






