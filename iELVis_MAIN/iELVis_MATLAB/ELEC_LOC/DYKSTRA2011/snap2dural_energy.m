function [coord_snapped,pairs]=snap2dural_energy(coord,surf,pairs)
% snap2dural_energy - uses energy minimization to snap electrodes to dural
%                     surface
%
% usage:
%   >> coord_snapped = snap2dural_energy(coord,surf);
%   >> [coord_snapped,pairs] = snap2dural_energy(coord,surf,pairs);
%
% inputs:
%   coord    Nx3 matrix of original RAS coordinates of electrodes
%   surf     Surface structure (single hemisphere) with .vert and .tri
%            fields.  Only pass in the surface of the appropriate
%            hemisphere (e.g. cortex.lh or cortex.rh)
%   pairs    Mx2 matrix of M pairs of electrode indices for which the deformation
%            energy is computed.  If not supplied, the 4 nearest neighbors
%            of each electrode are used.
%
% outputs:
%   coord_snapped  Nx3 matrix of coordinates snapped to dural surface
%   pairs          Mx2 matrix of the electrode pairs used in the
%                  deformation energy
%
% algorithm details: This function performs a constrained optimization
% using an energy function which accounts for the distance the electrodes
% move away from their starting point, and the change in distance of the
% electrodes from their neighbors.  The constraint that the electrodes need
% to lie on the dural surface is satisfied.


% compute pairs of neighbors
if(nargin<3 || isempty(pairs))
    pairs=knn_pairs(coord,4);
end

% set starting coordinates
coord0=coord;

% anonymous function handles
efun=@(coord_snapped)energy_electrodesnap(coord_snapped,coord,pairs);
cfun=@(coord_snapped)surface_constraint(coord_snapped,surf);

% options
options=optimset('Algorithm','active-set',...
                 'Display','iter',...
                 'MaxIter',50,...
                 'MaxFunEvals',Inf,...
                 'UseParallel','always',...
                 'GradObj','off',...
                 'TypicalX',reshape1d(coord),...
                 'DiffMaxChange',2,...
                 'DiffMinChange',0.3,...
                 'TolFun',0.3,...
                 'TolCon',0.01*size(coord0,1),...
                 'TolX',0.5,...
                 'Diagnostics','off',...
                 'RelLineSrchBnd',1,...
                 'PlotFcns',{@optimplotfval,@optimplotstepsize,@optimplotfirstorderopt,@plotCoordFun});
             
% history
history.coords=[];
last_handle=[];
             
% run minimization
coord_snapped = fmincon(efun,coord0,[],[],[],[],[],[],cfun,options);

% energy function
    function energy=energy_electrodesnap(coord,coord_orig,pairs)

        % energy from change in electrode positions
        energy_eshift=sum((coord-coord_orig).^2,2);
        
        % energy from deformation of grid shape
        dist=sqrt(sum((coord(pairs(:,1),:)-coord(pairs(:,2),:)).^2,2));
        dist_orig=sqrt(sum((coord_orig(pairs(:,1),:)-coord_orig(pairs(:,2),:)).^2,2));
        energy_deform=(dist-dist_orig).^2;

        % total energy
        energy=mean(energy_eshift)+10*mean(energy_deform.^2);
    end

% constraint function
    function [c,ceq]=surface_constraint(coord,surf)
        ind=knnsearch(coord,surf.vert,1);
        dist=sqrt(sum((coord-surf.vert(ind,:)).^2,2));

        c=[];
        ceq=dist;
    end

% nearest neighbor function
    function pairs=knn_pairs(coord,k)
        % compute pairs
        if(nargin<2 || isempty(k))
            k=5;
        end
        knn_ind=knnsearch(coord,coord,k);
        firstofpair=repmat([1:size(coord,1)]',1,k);
        kdist=reshape(sqrt(sum((coord(knn_ind(:),:)-coord(firstofpair(:),:)).^2,2)),size(knn_ind));
        
        pairs_all=cat(3,knn_ind,firstofpair);
        pairs_all=permute(pairs_all,[3 1 2]);
        pairs_all=sort(reshape(pairs_all,2,[]),1)';
        pairs_all=unique(pairs_all,'rows');
        
        dist_pairs=sqrt(sum((coord(pairs_all(:,1),:)-coord(pairs_all(:,2),:)).^2,2));
        dist_pairs=round(dist_pairs*20)/20;
        dist_standard=mode(dist_pairs);
        
        pairs=[];
        for i=1:size(knn_ind,1)
            numvalid=sum(kdist(i,:)<(dist_standard*1.25));
            if(numvalid==0)
                numvalid=sum(kdist(i,:)<(kdist(i,1)*1.25));
            end
            pairs=cat(1,pairs,[i*ones(numvalid,1),knn_ind(i,1:numvalid)']);
        end
        pairs=sort(pairs,2);
        pairs=unique(pairs,'rows');
    end

% nearest neighbor function old
    function pairs=knn_pairs_old(coord,k)
        % compute pairs
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

% plotting function
    function stop = plotCoordFun(coord,optimValues,state)

        switch state
            case 'init'
                hold on; axis vis3d; view([3 1 -0.5]);
                last_handle=scatter3(coord(:,1),coord(:,2),coord(:,3),'marker','o','markerfacecolor','r','markeredgecolor','r','sizedata',15);
                history.coords=coord;
            case 'iter'
                set(last_handle,'marker','.','sizedata',10,'markerfacecolor','b','markeredgecolor','b');
                last_handle=scatter3(coord(:,1),coord(:,2),coord(:,3),10,'marker','o','markerfacecolor','r','markeredgecolor','r','sizedata',15);
                plot3([coord(:,1) history.coords(:,1,end)]',[coord(:,2) history.coords(:,2,end)]',[coord(:,3) history.coords(:,3,end)]','k');
                history.coords=cat(3,history.coords,coord);
            case 'done'
                hold off;
            otherwise
        end
        stop = false;
    end

end