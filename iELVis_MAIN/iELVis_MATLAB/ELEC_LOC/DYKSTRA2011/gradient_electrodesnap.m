function dcoord=gradient_electrodesnap(coord,coord0,pairs)

dshift=2*(coord-coord0);

ddeform=zeros(size(coord));
for i=1:size(coord,1)
    [rind,cind]=find(pairs==i);
    cdiff=coord(pairs(rind,1),:)-coord(pairs(rind,2),:);
    c0diff=coord0(pairs(rind,1),:)-coord0(pairs(rind,2),:);
    energy=sum(cdiff.^2,2)-sum(c0diff.^2,2);
    ddeform(i,:)=4*sum(2*repmat(energy,1,3).*cdiff,1);
end

dcoord=20*dshift./size(coord,1)+ddeform./size(pairs,1);