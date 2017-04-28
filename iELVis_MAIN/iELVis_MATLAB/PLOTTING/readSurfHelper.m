function cort=readSurfHelper(surf_fname)
%function cort=readSurfHelper(surf_fname)
%
% Same as freesurfer's read_surf.m but adds a 1 to cort.tri if needed
%
% -David Groppe

[cort.vert cort.tri]=read_surf(surf_fname);
if min(min(cort.tri))<1
    cort.tri=cort.tri+1; %sometimes this is needed sometimes not. no comprendo. DG ??
end