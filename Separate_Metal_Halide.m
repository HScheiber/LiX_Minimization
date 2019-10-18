function [Metal,Halide] = Separate_Metal_Halide(Salt)
% Preprocessing
[match,~] = regexp(Salt,{'Li' 'Na' 'K' 'Rb' 'Cs' 'F' 'Cl' 'Br' 'I'},'match','tokens');
Matched_Metals = match(1:5);
Matched_Halides = match(6:9);
matches_metals = Matched_Metals(~cellfun('isempty',Matched_Metals));
matches_halides = Matched_Halides(~cellfun('isempty',Matched_Halides));
if ~isempty(matches_metals)
    Metal = matches_metals{1}{1};
else
    Metal = '';
end

if ~isempty(matches_halides)
    Halide = matches_halides{1}{1};
else
    Halide = '';
end
end