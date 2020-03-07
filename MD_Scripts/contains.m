function output = contains(str,pattern,extra,extra2)
	if strcmpi(extra,'ignorecase') && extra2
		match = regexpi(str,pattern,'ONCE');
	else
		match = regexp(str,pattern,'ONCE');
	end
    if isempty(match)
        output = false;
    else
        output = true;
    end
	
end