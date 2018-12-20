function [showElecCoords, showElecNames, showElecIds]=getElecPlottingInfoBids(surfType,coordType,bidsDir,fsSub,side,ignoreDepthElec,elecNames,bidsSes)
%function [showElecCoords, showElecNames, showElecIds]=getElecPlottingInfoBids(surfType,coordType,bidsDir,fsSub,side,ignoreDepthElec,elecNames,bidsSes)
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
        coordType='inf';
        warning('Using *.INF electrode coordinates as this is required for inflated brain.');
    end
end
% Possible coord types from plotPialSurf: 'LEPTO','POSTIMPLANT','PIAL', or 'INF'
% Corresponding iEEG-BIDS coord types: 'lepto', 'postimplant', 'pial',
% 'inf'
coordType=lower(coordType);
%/Users/davidgroppe/Desktop/HandMotor/sub-PT001/ieeg/sub-PT001_ses-01_space-inf_electrodes.tsv
coordFname=fullfile(bidsDir,['sub-' fsSub],'ieeg',sprintf('sub-%s_ses-%.2d_space-%s_electrodes.tsv',fsSub,bidsSes,coordType));
fprintf('Taking electrode info from %s. Use cfg.eleccord=''n''; if not wanted.\n',coordFname);

%% Collect electrode coordinates
elecCoordCsv=csv2Cell(coordFname,9,0); %9=tab
nElecTotal=size(elecCoordCsv,1)-1;
RAS_coor=zeros(nElecTotal,3);
coordHdrs={'x','y','z'};
for csvLoopB=1:3,
    colId=findStrInCell(coordHdrs{csvLoopB},elecCoordCsv(1,:),1);
    for csvLoopA=1:nElecTotal,
        RAS_coor(csvLoopA,csvLoopB)=str2double(elecCoordCsv{csvLoopA+1,colId});
    end
end

%% Collect electrode name, type, & hemisphere
elecInfo=cell(nElecTotal,3); % electrode names, hemisphere, and type
nameId=findStrInCell('name',elecCoordCsv(1,:),1);
typeId=findStrInCell('type',elecCoordCsv(1,:),1);  % grid, depth, strip
hemId=findStrInCell('hemisphere',elecCoordCsv(1,:),1); % L,R
for csvLoopA=1:nElecTotal,
    elecInfo{csvLoopA,1}=elecCoordCsv{csvLoopA+1,nameId}; % Name
    switch elecCoordCsv{csvLoopA+1,typeId}
        case 'grid'
            elecInfo{csvLoopA,2}='G';
        case 'strip'
            elecInfo{csvLoopA,2}='S';
        case 'depth'
            elecInfo{csvLoopA,2}='D';
        otherwise
            error('Unrecognized value of electrode "type": %s',elecCoordCsv{csvLoopA+1,typeId});
    end
    elecInfo{csvLoopA,3}=elecCoordCsv{csvLoopA+1,hemId}; % Hemisphere
end

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