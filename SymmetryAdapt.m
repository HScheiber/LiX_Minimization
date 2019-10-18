function [b,c] = SymmetryAdapt(a_in,b_in,c_in,Structure)
switch Structure
    case 'BetaBeO'
        %c = (1/sqrt(3)).*a_in;
        c = c_in;
        b = a_in;
        
    case 'CsCl'
        b = a_in;
        c = a_in;

    case 'FiveFive'
        %b = (3/2).*a_in; % Theoretical b/a = 3/2
        %c = (sqrt(3)/2).*a_in; % Theoretical c/a = cosd(30)
        b = b_in;
        c = c_in;
        
    case 'NiAs'
        %x = 1.391666666667;
        c = c_in; % x.*a_in; 1.7316 LiF LSDA, 1.7275 LiF PBE
        b = a_in;
    case 'Rocksalt'
        b = a_in;
        c = a_in;
        
    case 'Sphalerite'
        b = a_in;
        c = a_in;
        
    case 'Wurtzite'
        %c = sqrt(8/3).*a_in; % Perfect Wurtzite c/a = sqrt(8/3);
        c = c_in;
        b = a_in;
end

end