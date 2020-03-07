function output = contains(str,pattern,varargin)
	if nargin > 2
		if strcmpi(varargin{1},'ignorecase') && varargin{2}
			match = regexpi(str,pattern,'ONCE');
		else
			match = regexp(str,pattern,'ONCE');
		end
	else
		match = regexp(str,pattern,'ONCE');
	end
    if isempty(match)
        output = false;
    else
        output = true;
    end
	
end