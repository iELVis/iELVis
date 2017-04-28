%---------------------------------------------------
% FUNCTION: Y = get_loc_snap_mgh(electrodes,surface_path,side,surftype)
% INPUTS:   
%           electrodes = Nx3 matrix of original RAS coordinates of electrodes
%           
%           surface_path = path to smoothed pial surface files output from freesurfer (i.e. lh.pial.smoothed & rh.pial.smoothed)
%
%                    The lh and rh structures can be read in to MATLAB
%                    using the read_surf function (e.g.
%                    [cortex.lh.vert, cortex.lh.tri] =
%                    read_surf('lh.pial'), where lh.pial is the filename of
%                    the reconstructed left hemisphere surface.
%
%   	    side = 
%               'r' = right               or . . .
%               'l' = left                
% 
% OUTPUT:   New RAS electrode coordinates which have been "snapped" to the
%           dural surface
%---------------------------------------------------
% written by: Andrew Dykstra, MIT, 2009
% last updated 2010.03.31
%--------------------------------------------------- 

function Y = get_loc_snap_mgh(electrodes,surface_path,side,surftype)
%snaps "best guess" of electrodes to template


if nargin < 4
  surftype='pial-outer-smoothed';
end

%read smoothed surfaces into matlab
[cortex.lh.vert, cortex.lh.tri] = read_surf([surface_path '/lh.' surftype]);
[cortex.rh.vert, cortex.rh.tri] = read_surf([surface_path '/rh.' surftype]);

if side == 'l'
    vert_cart = cortex.lh.vert;
elseif side == 'r'
    vert_cart = cortex.rh.vert;
end

%Get spherical coordinates just in case
%[vert_sph(:,1), vert_sph(:,2), vert_sph(:,3)] = cart2sph(vert_cart(:,1),vert_cart(:,2),vert_cart(:,3));
%[elec_sph(:,1), elec_sph(:,2), elec_sph(:,3)] = cart2sph(elec_cart(:,1),elec_cart(:,2),elec_cart(:,3));

%loop over electrodes
for i = 1:max(size(electrodes))
    
    b_x = (vert_cart(:,1)-electrodes(i,1)).^2;
    b_y = (vert_cart(:,2)-electrodes(i,2)).^2;
    b_z = (vert_cart(:,3)-electrodes(i,3)).^2;
    
    temp = sqrt(b_x + b_y + b_z);
    [temp, index] = min(temp);
    clear temp*
    
    %Spherical Coordinates, assign electrode location to closest vertex
        %th_new(i) = vert_sph(index,1);
        %phi_new(i) = vert_sph(index,2);
        %r_new(i) = vert_sph(index,3)+1;
        
    %Cartesian Coordiantes, assign electrode location to closest vertex
        x_new(i) = vert_cart(index,1);
        y_new(i) = vert_cart(index,2);
        z_new(i) = vert_cart(index,3);
end
Y = [x_new' y_new' z_new'];

%if distance metric was computed in spherical coordiantes, convert back to
%cartesion
%[Y(:,1), Y(:,2), Y(:,3)] = sph2cart(th_new,phi_new,r_new);
