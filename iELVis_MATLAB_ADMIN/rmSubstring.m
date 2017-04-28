function outStr=rmSubstring(inStr,subStr,caseSensitive)
%function outStr=rmSubstring(inStr,subStr,caseSensitive)
%
% This function removes the first occurrence of a substring from a string.
%
% Inputs:
%   inStr - Input string
%   subStr - Substring to be removed
%
% Optional Input:
%   caseSensitive - [1 or 0] If 0, letter case will be ignored when finding
%                   a match of the substring {default: 0}.
%
% See also: rm_char.m
% David M. Groppe

if nargin<3,
    caseSensitive=0; 
end

if universalYes(caseSensitive)
    start_id=findstr(subStr,inStr);
else
    start_id=findstr(lower(subStr),lower(inStr));
end
if isempty(start_id)
    outStr=inStr;
else
    bad_ids=start_id:start_id+length(subStr)-1;
    outStr=inStr(setdiff(1:length(inStr),bad_ids));
end

