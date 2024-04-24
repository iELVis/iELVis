function [showElecColors, elecCbarMin, elecCbarMax, elecCmapName]=elec2rgb(elecColors,showElecIds,elecColorScale,elecCbar,elecCmapName)
%function [showElecColors, elecCbarMin, elecCbarMax, elecCmapName]=elec2rgb(elecColors,showElecIds,elecColorScale,elecCbar,elecCmapName)
%
% This function derives the RGB color values of electrodes. It gets called
% from plotElecs.m
%
% Inputs:
%  elecColors    -Empty | r | vector of elec values | nx3 matrix of rgb
%                 values (range 0-1); Electrode colors. Empty=all black. 'r'=all red.
%                 elec values are converted to rgb values; nx3 matrix fully
%                 specifies rgb values
%  showElecIds   -Indices of electrodes to show. May be a subset of
%                 elecColors
%  elecColorScale-see plotPialSurf.m documentation
%  elecCbar      -'y' or 'n'; if 'y' a colorbar will be plot for electrodes.
%
% Outputs:
% showElecColors-mx3 matrix of rgb values of electrodes that will be
%                visualized (might be a subset of elecColors)
% elecCbarMin   -max value of visualized electrodes
% elecCbarMax   -min value of visualized electrodes
% elecCmapName  -the colormap used for the electrodes

nShowElec=length(showElecIds);
elecCbarMin=[];
elecCbarMax=[];
%elecCmapName=[];
if isempty(elecColors)
    % No color/value info given, make electrodes black
    showElecColors = zeros(nShowElec,3);
elseif ischar(elecColors) && strcmp(elecColors,'r')
    % Make electrodes red
    showElecColors = zeros(nShowElec,3);
    showElecColors(:,1) = 1;
elseif isvector(elecColors) && size(elecColors,2)~=3
    % electrode values passed, automatically convert to colors
    % we need the second condition in case user wants to plot a single
    % electrode and passes an rgb vector to specify the color
    
    if isnumeric(elecColorScale)
        type='minmax';
        elecCbarMin=elecColorScale(1);
        elecCbarMax=elecColorScale(2);
    else
        type=elecColorScale;
    end
    if isempty(elecCmapName)
        if verLessThan('matlab','8.0.1')
            elecCmapName='jet';
        else
            elecCmapName='parula';
        end
    end
    if isempty(elecCbarMin),
        % make electrode colormap
        [showElecColors, elecLimits, elecCmapName]=vals2Colormap(elecColors(showElecIds),type,elecCmapName);
        elecCbarMin=elecLimits(1);
        elecCbarMax=elecLimits(2);
    else
        [showElecColors, elecLimits, elecCmapName]=vals2Colormap(elecColors(showElecIds),type,elecCmapName,[elecCbarMin elecCbarMax]);
        elecCbarMin=elecLimits(1);
        elecCbarMax=elecLimits(2);
    end
else
    showElecColors=elecColors(showElecIds,:);
%     if ~universalNo(elecCbar),
        if isnumeric(elecColorScale) && isvector(elecColorScale) && length(elecColorScale)==2
            elecCbarMin=min(elecColorScale);
            elecCbarMax=max(elecColorScale);
        else
            error('When cfg.elecColors is a matrix of RGB values, elecColorScale needs to specify the min and max of the colorscale.');
        end
%     else
%        elecCbarMin=[];
%        elecCbarMax=[];
%     end
end