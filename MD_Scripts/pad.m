function stringReturn = pad(stringPassed,totalChars,charPosition,fillChar)
if nargin == 0
    stringReturn = '';
    return;
end
if nargin<4
    fillChar = ' ';
    if nargin<3
        charPosition='right';
        if nargin<2
            stringReturn = stringPassed;
            return
        end
    end
end
if length(stringPassed)>=totalChars
    stringReturn = stringPassed;
    return
end
if size(fillChar,1) ~= 1 || size(fillChar,2) ~=1
    warning('The fill char pass is too large using space instead');     %#ok<WNTAG>
    fillChar = ' ';
end
% Go through from the current length to the desired length the required len
stringReturn = stringPassed;
for i=length(stringPassed)+1:totalChars
    if strcmp(charPosition,'left')
        stringReturn = [fillChar,stringReturn];     %#ok<AGROW>
    elseif strcmp(charPosition,'right')
        stringReturn = [stringReturn,fillChar];     %#ok<AGROW>
    end
end
end