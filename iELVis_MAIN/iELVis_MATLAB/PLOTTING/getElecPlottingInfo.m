function [showElecCoords, showElecNames, showElecIds]=getElecPlottingInfo(surfType,coordType,fsDir,fsSub,side,ignoreDepthElec,elecNames)
%function [showElecCoords, showElecNames, showElecIds]=getElecPlottingInfo(surfType,coordType,fsDir,fsSub,side,ignoreDepthElec,elecNames)
%
% This function reads electrode coordinates and names from the elec_recon
% subfolder of the subject's freesurfer folder. It should only be called by
% plotElecs.m
% 
% Inputs:
%  see plotPialSurf.m documentation
%
% Outputs:
%  showElecCoords-mx3 matrix of electrode coordinates
%  showElecNames -mD cell array of electrode names
%  showElecIds   -vector of indices of all electrodes indicating which
%                 electrodes are actually visualized

if strcmpi(surfType,'inflated')
    if ~strcmp(coordType,'INF'),
        coordType='INF';
        warning('Using *.INF electrode coordinates as this is required for inflated brain.');
    end
end
coordFname=fullfile(fsDir,fsSub,'elec_recon',[fsSub '.' coordType]);
fprintf('Taking electrode coordinates from %s.%s and %s.electrodeNames in elec_recon folder. Use cfg.eleccord=''n''; if not wanted.\n', ...
    fsSub,coordType,fsSub);
elecCoordCsv=csv2Cell(coordFname,' ',2);
nElecTotal=size(elecCoordCsv,1);
RAS_coor=zeros(nElecTotal,3);
for csvLoopA=1:nElecTotal,
    for csvLoopB=1:3,
        RAS_coor(csvLoopA,csvLoopB)=str2double(elecCoordCsv{csvLoopA,csvLoopB});
    end
end

elecInfoFname=fullfile(fsDir,fsSub,'elec_recon',[fsSub '.electrodeNames']);
elecInfo=csv2Cell(elecInfoFname,' ',2);

% Remove elecs in opposite hemisphere
if side=='l'
    showElecIds=find(cellfun(@(x) strcmp(x,'L'),elecInfo(:,3)));
else
    showElecIds=find(cellfun(@(x) strcmp(x,'R'),elecInfo(:,3)));
end

% Remove depth elecs if requested
if universalYes(ignoreDepthElec),
    fprintf('...not plotting depth electrodes (if any exist)\n');
    depthElecIds=find(cellfun(@(x) strcmp(x,'D'),elecInfo(:,2)));
    showElecIds=setdiff(showElecIds,depthElecIds);
end

allElecNames=elecInfo(:,1);
% If elecNames specified, remove any electrodes not in elecNames
rmIds=[];
if ~isempty(elecNames)
    for a=1:length(allElecNames),
        if isempty(findStrInCell(allElecNames{a},elecNames))
           rmIds=[rmIds a]; 
        end
    end
    showElecIds=setdiff(showElecIds,rmIds);
end
showElecCoords=RAS_coor(showElecIds,:);
showElecNames=allElecNames(showElecIds);