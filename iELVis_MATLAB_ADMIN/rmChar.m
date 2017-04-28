function out_str=rmChar(in_str,subchar)
%function out_str=rmChar(in_str,subchar)
%
% This function removes ALL occurrences of a character from a string
%
% See also: rm_substring.m
% David M. Groppe

bad_ids=find(subchar==in_str);
out_str=in_str(setdiff(1:length(in_str),bad_ids));
