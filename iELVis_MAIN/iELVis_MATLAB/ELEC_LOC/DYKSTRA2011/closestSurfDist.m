function dist=closestSurfDist(coord,surf)
ind=knnsearch(coord,surf.vert,1);
dist=sqrt(sum((coord-surf.vert(ind,:)).^2,2));
