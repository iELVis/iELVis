function yesno=universalYes(char_or_num)
%function yesno=universalYes(char_or_num)
%
% Returns 1 if char_or_num is 1, "on", "y", or "yes" and 0 otherwise.
% Useful for interpreting function arguments
%

if isempty(char_or_num)
    yesno=0;
elseif isnumeric(char_or_num) && (char_or_num==1)
    yesno=1;
elseif strcmpi(char_or_num,'yes') || ...
        strcmpi(char_or_num,'on') || strcmpi(char_or_num,'y')
    yesno=1;
else
    yesno=0;
end