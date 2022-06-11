function [showElecCoords, showElecNames, h_elec, elecCbarMin, elecCbarMax, elecCmapName]=plotElecs(elecCoord,surfType,fsDir,fsSub,side,ignoreDepthElec,pullOut,elecColors,elecColorScale,elecShape,elecSize,showLabels,clickElec,elecAssign,edgeBlack,edgeColors,elecNames,elecCbar,bidsDir,bidsSes)
% function [showElecCoords, showElecNames, h_elec, elecCbarMin, elecCbarMax, elecCmapName]=plotElecs(elecCoord,surfType,fsDir,fsSub,side,ignoreDepthElec,pullOut,elecColors,elecColorScale,elecShape,elecSize,showLabels,clickElec,elecAssign,edgeBlack,edgeColors,elecNames,elecCbar,bidsDir,bidsSes)
% This function plots electrodes. It should only be called by plotPialSurf.m
%
% Inputs:
%   see plotPialSurf.m documentation
%
% Outputs:
%  showElecCoords - nx3 matrix; The coordinates of the visualized electrodes
%  showElecNames  - cell array; The names of the visualized electrodes
%  h_elec         - cell array; handles of visualized electrodes
%  elecCbarMin    - scalar; min value of electrode colorbar (empty if not
%                   used)
%  elecCbarMax    - scalar; max value of electrode colorbar (empty if not
%                   used)
%

% Copy original value of elecColors to be able to tell if colored
% electrodes were desired even if no electrodes get shown
elecColorsOrig=elecColors;
elecCmapName=[];

% Get electrode coordinates
if isnumeric(elecCoord)
    % Electrode coordinates passed as argument
    if size(elecCoord,2)==4
        display('...Electrode input is matrix with coordinates.');
        showElecCoords=elecCoord;
        if side=='l',
            showElecIds=find(elecCoord(:,4));
        else
            showElecIds=find(~elecCoord(:,4));
        end
        
        % Remove electrodes over hemisphere not being shown
        showElecNames=elecNames(showElecIds);
        showElecCoords=showElecCoords(showElecIds,:);
        
        if strcmpi(surfType,'inflated')
            % Plotting on inflated brain. Convert electrode coordinates to
            % inflated coordinates
            cfg_pvox2inf=[];
            cfg_pvox2inf.fsurfSubDir=fsDir;
            cfg_pvox2inf.elecCoord=showElecCoords;
            cfg_pvox2inf.elecNames=showElecNames;
            showElecCoords=pial2InfBrain(fsSub,cfg_pvox2inf);
        end
        showElecCoords=showElecCoords(:,1:3); % Remove 4th column that indicates hemisphere each elec belongs to
    else
        error('elecCoord is numeric but doesn''t have 3 coordinates + binary hemisphere column. It needs to be an nx4 matrix where n=# of electrodes.');
    end
else
    if isempty(bidsDir)
        % Load electrode coordinates from subject's Freesurfer folder
        [showElecCoords, showElecNames, showElecIds]=getElecPlottingInfo(surfType,elecCoord,fsDir,fsSub,side,ignoreDepthElec,elecNames);
    else
        % Load electrode coordinates from subject's iEEG-BIDS directory
        [showElecCoords, showElecNames, showElecIds]=getElecPlottingInfoBids(surfType,elecCoord,bidsDir,fsSub,side,ignoreDepthElec,elecNames,bidsSes);
    end
    if ~isempty(elecColors) && isnumeric(elecColors),
        % If elec colors/values have been specified, expand them so that
        % they correspond to showElecIds, since colors/values might only be
        % specified at a subset of electrodes
        n_col=size(elecColors,2);
        tempElecColors=zeros(max(showElecIds),n_col);
        for a=1:length(elecNames),
            temp_id=findStrInCell(elecNames{a},showElecNames);
            if ~isempty(temp_id),
                tempElecColors(showElecIds(temp_id),:)=elecColors(a,:);
            end
        end
        elecColors=tempElecColors;
    end
end

% Pull electrodes out from the brain towards the viewer
nShowElec=size(showElecCoords,1);
h_elec=cell(nShowElec); %preallocate memory
if pullOut,
    fprintf('...pulling out electrodes by factor %f. cfg.pullOut=0 if not wanted.\n',pullOut);
    v=axis;
    campos=get(gca,'cameraposition');
    %camtarg=get(gca,'cameratarget'); ?? Should we check that this is set to 0?
    err=repmat(campos,nShowElec,1)-showElecCoords;
    nrmd=err./repmat(sqrt(sum(err.^2,2)),1,3);
    showElecCoords=showElecCoords+nrmd*pullOut;
end

% Figure out electrode colors
elecCbarMin=[];
elecCbarMax=[];
[showElecColors, elecCbarMin, elecCbarMax, elecCmapName]=elec2rgb(elecColors,showElecIds,elecColorScale,elecCbar);
if ~isempty(elecColorsOrig) && isempty(showElecColors),
   error('ERROR: You tried to plot colored electrodes, but none of them are visible.'); 
end

% Prepare variables if electrodes are to be drawn as spheres
if strcmpi(elecShape,'sphere')
    elecSphere=1;
    [sphX, sphY, sphZ]=sphere(20);
    Zdim=size(sphZ);
    scale_sph=elecSize;
    sphX=sphX*scale_sph;
    sphY=sphY*scale_sph;
    sphZ=sphZ*scale_sph;
    sph_colors=zeros(nShowElec,3);
    sph_ct=0;
else
    elecSphere=0;
end

% Plot Electrodes
for j = 1:nShowElec
    if elecSphere
        sph_ct=sph_ct+1;
        h_elec{sph_ct}=surf(sphX+showElecCoords(j,1),sphY+showElecCoords(j,2),sphZ+showElecCoords(j,3),zeros(Zdim));
        sph_colors(sph_ct,:)=showElecColors(j,:);
    else
        h_elec{j}=plot3(showElecCoords(j,1),showElecCoords(j,2),showElecCoords(j,3),'o','Color',showElecColors(j,:),'MarkerFaceColor', showElecColors(j,:),'MarkerSize',elecSize);
        
        if (iscell(edgeBlack) && max(size(edgeBlack))==1) && universalYes(edgeBlack)...
                || (ischar(edgeBlack) && universalYes(edgeBlack))
            set(h_elec{j},'MarkerEdgeColor','k');
        elseif any(strcmp(showElecNames{j},edgeBlack))
            set(h_elec{j},'MarkerEdgeColor','k');
        else
            set(h_elec{j},'MarkerEdgeColor',showElecColors(j,:));
        end
        
        if ~isempty(edgeColors)
            [ind1, ~] = match_str(edgeColors(:,1),showElecNames{j});
            if ~isempty(ind1)
               set(h_elec{j},'MarkerEdgeColor',edgeColors{ind1,2});
            end          
        end
        
        
        
        if showLabels
            add_name(showElecCoords(j,:),showElecNames{j},showElecNames,elecSize,showElecColors(j,:))
        end
    end
    hold all
    if universalYes(clickElec),
        if isempty(elecAssign)
            if elecSphere,
                set(h_elec{sph_ct},'userdata',showElecNames{j});
            else
                set(h_elec{j},'userdata',showElecNames{j});
            end
        else
            if elecSphere,
                set(h_elec{sph_ct},'userdata',showElecNames{j});
            else
                set(h_elec{j},'userdata',[showElecNames{j} ' ' elecAssign{j,2}]);
            end
        end
        % This click_text code should put the text out towards the
        % viewer (so it doesn't get stuck in the brain)
        % Note: pop_fact=5 in the below code might be too far for lateral surfaces
        bdfcn=['Cp = get(gca,''CurrentPoint''); ' ...
            'Cp=Cp(1,1:3);', ...
            'v=axis;', ...
            'campos=get(gca,''cameraposition'');', ...
            'df=Cp-campos;', ...
            'nrmd=df/sqrt(sum(df.^2));', ...
            'pop_fact=5;', ...
            'eval(sprintf(''Cp=Cp-%d*nrmd;'',pop_fact));', ...
            'dat=get(gcbo,''userdata'');', ...
            'ht=text(Cp(1),Cp(2),Cp(3),sprintf(''%s'',dat));', ...
            'set(ht,''backgroundColor'',''w'',''horizontalalignment'',''center'',''verticalalignment'',''middle'',''buttondownfcn'',''delete(gcbo);'');'];
        if elecSphere,
            set(h_elec{sph_ct},'buttondownfcn',bdfcn);
        else
            set(h_elec{j},'buttondownfcn',bdfcn);
        end
    end
end
%NOTE:
% x dimension is lateral/medial  (+=lateral)
% y dimension is ant/posterior (+=anterior)
% z dimension is superior/inferior (+=superior)

if elecSphere,
    shading interp; lighting gouraud; material dull;
    %for some reason the shading command resets of the colors of all the
    %spheres, thus the need for this silly loop.  There is surely a more
    %elegant way to deal with this.
    for a=1:sph_ct,
        set(h_elec{a},'facecolor',sph_colors(a,:));
    end
end



%% subfunction add_name
function add_name(xyz,label,all_labels,markersize,rgb)
% Adds an electrode's name next to its location
% right now rgb argument is ignored, fix in the future so that electrode
% names stand out from background? ??

if keyElec(label,all_labels),
    h_t=text(xyz(1),xyz(2),xyz(3),label);
    set(h_t,'color','k','fontweight','bold','fontsize',markersize+4);
    %     if isequal(rgb,[0 0 0]),
    %         set(h_t,'color','w','fontweight','bold','fontsize',markersize+2);
    %     else
    %         set(h_t,'color','k','fontweight','bold','fontsize',markersize+2);
    %     end
end
