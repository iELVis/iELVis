function energy=deformation_energy(coord,coord_orig,pairs)

dist=sqrt(sum((coord(pairs(:,1),:)-coord(pairs(:,2),:)).^2,2));
dist_orig=sqrt(sum((coord_orig(pairs(:,1),:)-coord_orig(pairs(:,2),:)).^2,2));

energy=(dist-dist_orig).^2;