function [energy denergy]=energy_electrodesnap(coord,coord_orig,pairs)

% energy_surf=closestSurfDist(coord,surf);
energy_eshift=sum((coord-coord_orig).^2,2);
% energy_deform=deformation_energy_full(coord,coord_orig);
energy_deform=deformation_energy(coord,coord_orig,pairs);

% energy=2*mean(energy_surf)+mean(energy_eshift)+0.5*mean(energy_deform);
energy=mean(energy_eshift)+mean(energy_deform.^2);
% fprintf('shift: %0.4f, deform: %0.4f, energy: %0.4f\n',mean(energy_eshift),mean(energy_deform),mean(energy));
% energy=mean(energy_deform);

% denergy=gradient_electrodesnap(coord,coord_orig,pairs);
denergy=[];