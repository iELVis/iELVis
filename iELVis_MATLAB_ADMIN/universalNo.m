function yesno=universalNo(char_or_num)
%function yesno=universalNo(char_or_num)
%
% Returns 1 if char_or_num is empty, 0, "off", "n", or "no" and 0 otherwise.
% Useful for interpreting function arguments
%

if ~isempty(char_or_num)
    if isnumeric(char_or_num) && isequal(char_or_num,0)
        yesno=1;
    elseif strcmpi(char_or_num,'no') || ...
            strcmpi(char_or_num,'off') || strcmpi(char_or_num,'n')
        yesno=1;
    else
        yesno=0;
    end
else
    yesno=1;
end