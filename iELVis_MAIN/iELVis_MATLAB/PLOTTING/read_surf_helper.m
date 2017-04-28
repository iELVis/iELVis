function cort=read_surf_helper(surf_fname)
%function cort=read_surf_helper(surf_fname)
%
% David Groppe

[cort.vert cort.tri]=read_surf(surf_fname);
if min(min(cort.tri))<1
    cort.tri=cort.tri+1; %sometimes this is needed sometimes not. no comprendo. DG ??
end