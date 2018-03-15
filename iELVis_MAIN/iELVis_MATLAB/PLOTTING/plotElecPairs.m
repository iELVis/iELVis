function plotElecPairs(elecPairs,lineWidth,side,elecNames,elecCoord,elecSize)
% plotElecPairs(elecPairs,lineWidth,side,elecNames,elecCoord,elecSize)
%
% Plots lines between electrodes. Useful for plotting bipolar data or
% physical connections between electrodes. This function is only used by
% plotPialSurf.m
%
% Inputs:
%  see plotPialSurf.m documentation

% Plot lines joining electrodes ?? update this
if ~isempty(elecPairs),
    % Remove all electrode pairs not on this hemisphere
    nPairs=length(elecPairs);
    usePairs=zeros(nPairs,1);
    for pairLoop=1:nPairs,
        if strcmpi(elecPairs{pairLoop,4},side)
            usePairs(pairLoop)=1;
        end
    end
    elecPairs=elecPairs(find(usePairs),:);
    clear usePairs nPairs
    
    if isempty(lineWidth)
        lineWidth=elecSize/3;
    end
    if size(elecPairs,2)<=5
        elecPairs(:,6) ={lineWidth};
    elseif size(elecPairs,2)>5
        % normalize lineWidth
        elecPairs(:,6) = cellfun(@rdivide,elecPairs(:,6), ...
            num2cell(repmat(max([elecPairs{:,6}]), [size(elecPairs,1) 1])),'UniformOutput',false);
        elecPairs(:,6) = cellfun(@times,elecPairs(:,6), ...
            num2cell(repmat(lineWidth, [size(elecPairs,1) 1])),'UniformOutput',false);
    end
    
    n_pairs=size(elecPairs,1);
    pair_ids=[0 0];
    for a=1:n_pairs,
        for b=1:2,
            [got_it, pair_ids(b)]=ismember(lower(elecPairs{a,b}),lower(elecNames));
            if ~got_it
                error('Channel %s is in electrode pairs but not in pialVox electrode names.',elecPairs{a,b});
            end
        end
        hl=plot3([elecCoord(pair_ids(1),1) elecCoord(pair_ids(2),1)], ...
            [elecCoord(pair_ids(1),2) elecCoord(pair_ids(2),2)], ...
            [elecCoord(pair_ids(1),3) elecCoord(pair_ids(2),3)],'-');
        if isnumeric(elecPairs{a,3})
            set(hl,'color',elecPairs{a,3},'lineWidth',lineWidth);
        else
            set(hl,'color',str2num(elecPairs{a,3}),'lineWidth',lineWidth);
        end
        
        if size(elecPairs,2)>4 && ~isempty(elecPairs{a,5})
            clickText3D(hl,[elecPairs{a,1} '-' elecPairs{a,2} ': ' elecPairs{a,5}],2);
        else
            clickText3D(hl,[elecPairs{a,1} '-' elecPairs{a,2}],2);
        end
    end
end