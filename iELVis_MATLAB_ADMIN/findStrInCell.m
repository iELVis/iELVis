function ids=findstrInCell(str,cell_o_strings,error_flag)
%function ids=findstrInCell(str,cell_o_strings,error_flag)
%
% Returns the indices of a target string in cell array of strings.
%
% Required Inputs:
%   str            - Target string
%   cell_o_strings - Cell array of strings
%
% Optional Input:
%   error_flag - [1|0] If 1, an error is thrown when str is not found in
%                cell_o_strings
%
% Output:
%   ids - Indices of elements of cell_o_strings that match str.  ids is
%         empty if there are no matches.
%
% Note, NOT case sensitive.
%
% Author: David Groppe

if nargin<3,
   error_flag='n'; 
end

ids=[];
n=length(cell_o_strings);

for a=1:n,
   if strcmpi(str,cell_o_strings{a})
      ids=[ids a]; 
   end
end

if universalYes(error_flag) && isempty(ids)
   error('String %s not found in cell array of strings.',str); 
end