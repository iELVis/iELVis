function pairs=knn_pairs(coord,k)
    if(nargin<2 || isempty(k))
        k=4;
    end
    knn_ind=knnsearch(coord,coord,k);
    pairs=cat(3,knn_ind,repmat([1:size(coord,1)]',1,k));
    pairs=permute(pairs,[3 1 2]);
    pairs=sort(reshape(pairs,2,[]),1)';
    pairs=unique(pairs,'rows');

    % compute threshold
    dist=pdist(coord,'euclidean');
    dist=sort(dist,'ascend');
    threshold=mean(dist(1:(size(coord,2)*3)));
    threshold=threshold*1.5;

    % take pairs less than threshold
    dist_pairs=sqrt(sum((coord(pairs(:,1),:)-coord(pairs(:,2),:)).^2,2));
    pairs=pairs(dist_pairs<=threshold,:);
end