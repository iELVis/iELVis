function [surf] = fs_load_subj(subj,hemi,surfname,nverts_only,subjdir)
%function [surf] = fs_load_subj(subj,hemi,[surfname],[nverts_only],[subjdir])
%
% Required input:
%  subj is a string specifying the subject name
%  hemi should be either 'lh' for left hemisphere or 'rh' for right hemi
%
% Optional parameters:
%  surfname - surface file to be loaded
%    {default: white}
%  nverts_only  - if 1, don't actually load surface, just get number of vertices
%    {default = 0}
%  subjdir - subjects directory (override SUBJECTS_DIR environment variable)
%    subjdir/subj should contain the freesurfer subject directory
%    {default = $SUBJECTS_DIR}
%
% Output:
%   surf is a structure containing:
%    nverts: number of vertices
%    nfaces: number of faces (triangles)
%    faces:  vertex numbers for each face (3 corners)
%    vertices: x,y,z coordinates for each vertex
%    nbrs:   vertex numbers of neighbors for each vertex
%
% created:        06/11/06 Don Hagler
% last modified:  02/24/10 Don Hagler
%
% see also: fs_read_surf()
%

if nargin < 2, help(mfilename); return; end;

surf = [];

if ~exist('surfname','var') || isempty(surfname)
  surfname = 'white';
end;

if ~exist('nverts_only','var') || isempty(nverts_only)
  nverts_only = 0;
end;

if ~ismember(hemi,{'lh','rh'})
  error('hemi must be lh or rh (is %s)',hemi);
end;

if ~exist('subjdir','var') || isempty(subjdir)
  subjdir = getenv('SUBJECTS_DIR');
  if isempty(subjdir)
    error('SUBJECTS_DIR not defined as an environment variable');
  end;
end;

surffile = sprintf('%s/%s/surf/%s.%s',subjdir,subj,hemi,surfname);
if ~exist(surffile,'file')
  error('surface file %s not found',surffile);
end

if nverts_only
  surf = fs_read_surf_nverts(surffile);
else
  surf = fs_read_surf(surffile);
  surf = fs_find_neighbors(surf);
end;

return;

