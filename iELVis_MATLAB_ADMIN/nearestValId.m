function nearestId=nearestValId(vec,scalarVal)
%function val=nearestValId(vec,scalarVal)
%
% Returns the value in the vector, vec, that is closest to the scalar,
% scalarVal.

df=abs(vec-scalarVal);
[~, nearestId]=min(df);