function [c,ceq]=surface_constraint(coord,surf)
c=[];
ceq=closestSurfDist(coord,surf);