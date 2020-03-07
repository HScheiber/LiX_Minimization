function output = contains(str,pattern)
    match = regexp(str,pattern,'ONCE');
    if isempty(match)
        output = false;
    else
        output = true;
    end
end