function [surf] = fs_read_surf(fname)
% fs_read_surf - read a freesurfer surface file
% 
% [surf] = fs_read_surf(fname)
% 
% Reads the vertex coordinates (mm) and face lists from a surface file
%
% surf is a structure containing:
%   nverts: number of vertices
%   nfaces: number of faces (triangles)
%   faces:  vertex numbers for each face (3 corners)
%   vertices: x,y,z coordinates for each vertex
%
% this is a modification of Darren Weber's freesurfer_read_surf
%   which was a modified version of freesurfer's read_surf
%
% see also: fs_read_trisurf, fs_find_neighbors, fs_calc_triarea
%
% created:        03/02/06 Don Hagler
% last modified:  03/31/10 Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%QUAD_FILE_MAGIC_NUMBER =  (-1 & 0x00ffffff) ;
%NEW_QUAD_FILE_MAGIC_NUMBER =  (-3 & 0x00ffffff) ;

TRIANGLE_FILE_MAGIC_NUMBER  =  16777214 ;
QUAD_FILE_MAGIC_NUMBER      =  16777215 ;

% open it as a big-endian file
fid = fopen(fname, 'rb', 'b') ;
if (fid < 0),
  error(sprintf('could not open surface file %s.',fname));
end

magic = fs_fread3(fid) ;

if (magic == QUAD_FILE_MAGIC_NUMBER),
  surf.nverts = fs_fread3(fid) ;
  surf.nfaces = fs_fread3(fid) ;
%  fprintf('%s: reading %d quad file vertices...',mfilename,surf.nverts); tic;
  surf.vertices = fread(fid, surf.nverts*3, 'int16') ./ 100 ; 
%  t=toc; fprintf('done (%0.2f sec)\n',t);
%  fprintf('%s: reading %d quad file faces (please wait)...\n',...
%    mfilename,surf.nfaces); tic;
  surf.faces = zeros(surf.nfaces,4);
  for iface = 1:surf.nfaces,
    for n=1:4,
      surf.faces(iface,n) = fs_fread3(fid) ;
    end
%    if(~rem(iface, 10000)), fprintf(' %7.0f',iface); end
%    if(~rem(iface,100000)), fprintf('\n'); end
  end
%  t=toc; fprintf('\ndone (%0.2f sec)\n',t);
elseif (magic == TRIANGLE_FILE_MAGIC_NUMBER),
%  fprintf('%s: reading triangle file...',mfilename); tic;

  tline = fgets(fid); % read creation date text line
  tline = fgets(fid); % read info text line

  surf.nverts = fread(fid, 1, 'int32') ; % number of vertices
  surf.nfaces = fread(fid, 1, 'int32') ; % number of faces

  % vertices are read in column format and reshaped below
  surf.vertices = fread(fid, surf.nverts*3, 'float32');

  % faces are read in column format and reshaped
  surf.faces = fread(fid, surf.nfaces*3, 'int32') ;
  surf.faces = reshape(surf.faces, 3, surf.nfaces)' ;
%  t=toc; fprintf('done (%0.2f sec)\n',t);
else
  error(sprintf('unknown magic number in surface file %s.',fname));
end

surf.vertices = reshape(surf.vertices, 3, surf.nverts)' ;
fclose(fid) ;

%fprintf('...adding 1 to face indices for matlab compatibility.\n\n');
surf.faces = surf.faces + 1;

return
