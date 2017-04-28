function grid_lines=derive_grid_lines(grid_stem,grid_dim,offset,verbLevel)
%function grid_lines=derive_grid_lines(grid_stem,grid_dim,offset,verbLevel)
%
% Required Inputs:
%  grid_stem - (string) The name of the grid (e.g., 'LGd')
%  grid_dim  - (two integer vector) The number of rows and columns in the
%              grid (e.g., [8 8]).
%
% Optional Inputs:
%  offset    - (integer) All electrode numbers will be incremented by this
%              integer. Useful for cut grids.
%  verbLevel - [1 or 0] If non-zero, the vectors of contact numbers that
%              are neighbors will be printed on the command line.
%
% Output:
%  grid_lines - (2D cell array) The first column is the grid_stem, the
%               second is a vector of neighbors in a row or column in the
%               grid.
%
% For use with get_bipolar_nbors.m
%
% Example:
% >grid_lines=derive_grid_lines('LGd',[8 8]);
%
% Author:
% David Groppe
% 2/2014
%

if nargin<3,
    offset=0;
end

if nargin<4,
    verbLevel=0;
end

grid_lines=[];
% Do rows
for a=1:grid_dim(1)
    grid_lines{a,1}=grid_stem;
    grid_lines{a,2}=[1:grid_dim(2)]+grid_dim(2)*(a-1)+offset;
    if universalYes(verbLevel)
        disp(grid_lines{a,2});
    end
end

%% Do columns
for a=1:grid_dim(2)
    grid_lines{grid_dim(1)+a,1}=grid_stem;
    grid_lines{grid_dim(1)+a,2}=a+[0:(grid_dim(1)-1)]*grid_dim(2)+offset;
    if universalYes(verbLevel)
        disp(grid_lines{grid_dim(1)+a,2});
    end
end
