function [vertex_coords, faces] = read_surf(fname)

%
% [vertex_coords, faces] = read_surf(fname)
% reads a the vertex coordinates and face lists from a surface file
% note that reading the faces from a quad file can take a very long
% time due to the goofy format that they are stored in. If the faces
% output variable is not specified, they will not be read so it
% should execute pretty quickly.
%


%fid = fopen(fname, 'r') ;
%nvertices = fscanf(fid, '%d', 1);
%all = fscanf(fid, '%d %f %f %f %f\n', [5, nvertices]) ;
%curv = all(5, :)' ;

% open it as a big-endian file


%QUAD_FILE_MAGIC_NUMBER =  (-1 & 0x00ffffff) ;
%NEW_QUAD_FILE_MAGIC_NUMBER =  (-3 & 0x00ffffff) ;

TRIANGLE_FILE_MAGIC_NUMBER =  16777214 ;
QUAD_FILE_MAGIC_NUMBER =  16777215 ;

fid = fopen(fname, 'rb', 'b') ;
if (fid < 0)
    str = sprintf('could not open curvature file %s.', fname) ;
    error(str) ;
   %oops(159) ;
end
magic = fread3(fid) ;

if (magic == QUAD_FILE_MAGIC_NUMBER)
    vnum = fread3(fid) ; %number of vertices
    fnum = fread3(fid) ; %number of faces
    vertex_coords = fread(fid, vnum*3, 'int16') ./ 100 ;
    if (nargout > 1)
         for i=1:fnum
             for n=1:4
                     faces(i,n) = fread3(fid) ;
            end
        end
    end
elseif (magic == TRIANGLE_FILE_MAGIC_NUMBER)
  fgets(fid) ;
    fgets(fid) ;
    vnum = fread(fid, 1, 'int32') ; %number of vertices
    fnum = fread(fid, 1, 'int32') ; %number of faces
    vertex_coords = fread(fid, vnum*3, 'float32') + 1 ;
    %The next line was added by Andrew Dykstra (HST-SHBT), Feb 12, 2009, in
    %order to give 3D coordinates for the vertices
    vertex_coords = reshape(vertex_coords, 3, vnum)';
    faces = fread(fid, fnum*3, 'int32') + 1 ;
    faces = reshape(faces, 3, fnum)' ;
end
fclose(fid) ;