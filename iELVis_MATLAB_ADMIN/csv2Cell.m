% csv2Cell() - Imports the contents of a character delimited spreadsheet 
%              text file as a 2D cell array, csv_matrix.  
%
% Usage:
%  >> csv_matrix=csv2Cell(csv_fname,delimiter,header_lines);
%
% Required Input:
%   csv_fname - The name of the text file to import.  If a path is not
%               included in the name, the current working directory will 
%               be used.
%
% Optional Inputs:
%   delimiter    - The character that marks cell boundaries in the text file.
%                  {default: ','}
%
%   header_lines - (integer) The first header_lines number of lines of the 
%                  csv file will not be imported.  Useful for ignoring 
%                  header rows. {default: 0}
%
%
% Output:
%   csv_matrix - 2D cell array of strings. csv_matrix{x,y} corresponds to
%                the cell in the xth row and yth column in the spreadsheet.
%
% Note, be sure not to use the delimiting character in the text file as 
% anything but a cell boundary marker.
%
% Example (comma delimited):
%  >> demo_dat=round(rand(5,3)*10);
%  >> fid=fopen('demo.csv','w');
%  >> for a=1:5, for b=1:3, fprintf(fid,'%d,',demo_dat(a,b)); end; fprintf(fid,'\n'); end;
%  >> fclose(fid);
%  >> csv_matrix=csv2Cell('demo.csv',',');
%
% Example (tab delimited):
%  >> demo_dat=round(rand(5,3)*10);
%  >> fid=fopen('demo.csv','w');
%  >> for a=1:5, for b=1:3, fprintf(fid,'%d\t',demo_dat(a,b)); end; fprintf(fid,'\n'); end;
%  >> fclose(fid);
%  >> csv_matrix=csv2Cell('demo.csv',9);
%
% Author: 
%  David Groppe
%  9/3/2011
%  Laboratory for Multimodal Human Brain Mapping
%  Feinstein Institute for Medical Research
%  Manhasset, New York

function csv_matrix=csv2Cell(csv_fname,delimiter,header_lines)

if nargin<1,
    help csv2Cell
end

if nargin<2,
    delimiter=',';
else
    if length(delimiter)>1
        error('Delimiter needs to be a single character.');
    elseif (delimiter<0) || (delimiter>127)
        error('Specified delimiter needs to be a character.');
    end
end

[fid, msg]=fopen(csv_fname,'r');
if fid==-1,
   error('Cannot open %s because: %s.\n',csv_fname,msg); 
end

if nargin<3,
   header_lines=0; 
end

for a=1:header_lines,
    tline = fgetl(fid);
end

csv_matrix=cell(1,1);
row_ct=1;
col_ct=1;
first_row_passed=0;
while 1
    tline = fgetl(fid);
    if ~ischar(tline), break, end

    while ~isempty(tline),
        [t, tline]=parse_by_char(tline,delimiter);
        csv_matrix{row_ct,col_ct}=t;
        col_ct=col_ct+1;
    end
    col_ct=1;

    row_ct=row_ct+1;
end
fclose(fid);


function [pre, post]=parse_by_char(str,delimiter)

if isspace(delimiter), char_ids=find(isstrprop(str,'wspace')); else, char_ids=find(str==delimiter); end % edited by Pierre M�gevand 2024/07/16 to deal with tabs as separators
if isempty(char_ids),
    pre=str;
    post=[];
else
    pre=str(1:char_ids(1)-1);
    post=str(char_ids(1)+1:end);
end