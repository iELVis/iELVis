function energy=deformation_energy_full(coord,coord_orig)

energy=abs(pdist(coord,'euclidean')-pdist(coord_orig,'euclidean'))./(pdist(coord_orig,'euclidean').^2);