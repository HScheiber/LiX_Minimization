function output = copy_atom_order(filename)
[~,~,ext] = fileparts(filename);

fid = fopen(filename,'rt');

if strcmp(ext,'.gro')
    Input_data = textscan(fid,'%5d%-5s%5s%5d%8.3f%8.3f%8.3f%8.4f%8.4f%8.4f\n',...
        'Delimiter','','whitespace','','headerlines',2);
    fclose(fid);

    N = length(Input_data{1,1})-1;
    Outcell = cell(1,N);
    for i=1:N  
        Outcell{i} = [Input_data{2}{i} ' 1'];
    end
elseif strcmp(ext,'.g96')
    curline = fgetl(fid);
    while ~strcmp(curline,'POSITION')
        curline = fgetl(fid);
    end
    Input_data = textscan(fid,'%6c%6c%6c%6c%15.9f%15.9f%15.9f\n',...
        'Delimiter','','whitespace','');
    fclose(fid);
    
    N = length(Input_data{3})-1;
    Outcell = cell(1,N);
    for i=1:N  
        Outcell{i} = [Input_data{3}(i,:) ' 1'];
    end
end

W = [Outcell',[repmat({newline},numel(Outcell)-1,1);{[]}]]';
output = [W{:}];
end