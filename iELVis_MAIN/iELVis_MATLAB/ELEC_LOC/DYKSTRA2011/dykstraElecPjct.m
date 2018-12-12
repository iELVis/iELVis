function dykstraElecPjct(sub,minimizeChange)
%function dykstraElecPjct(sub,minimizeChange)
%
% Corrects intracranial electrode locations for brain shift using the
% following method:
%  Dykstra, A.R., Chan, A.M., Quinn, B.T., Zepeda, R., Keller, C.J.,
%  Cormier, J., Madsen, J.R., Eskandar, E.N., Cash, S.S., 2011.
%  Individualized localization and cortical surface-based registration of
%  intracranial electrodes. NeuroImage 1-42.
%
% Required Input:
%  sub - Freesurfer subject name (e.g., 'TWH001')
%
% Optional Input:
%  minimizeChange - [0 | 1] If 0, subdural electrodes are simply projected
%        to the closest leptomeningeal surface vertex to correct for brain shift. If
%        non-zero, an interative optimization algorithm is used to project
%        electrodes to the leptomeningeal surface in a way tha minimizes their
%        change in spatial location and change in the Euclidean distance 
%        between electrode neighbors. Turning off the optimization may help 
%        if the optimization peforms poorly. {default: 1}
%
% Outputs:
%  The following files are created in the elec_recon subfolder of the
%  Freesufer subject folder:
%     *.POSTIMPLANT: The RAS coordinates of electrodes before any correction for 
%           postimplant brain shift
%     *.LEPTO: The leptomeningeal surface RAS coordinates of electrodes after 
%           correction for postimplant brain shift. Depth electrode 
%           coordinates are the same as in *.CT
%     *.LEPTOVOX: The leptomeningeal surface voxel coordinates of electrodes after 
%           correction for postimplant brain shift. Voxel coordinates 
%           are for brainmask.nii.gz file also in the elec_recon folder. 
%           Depth electrode coordinates are the same as in *.CT
%     *.INF: The pial RAS coordinates of electrodes on the inflated pial 
%           surface after correcting for brain shift. Depths have NaN coordinates.
%     *.PIAL: The pial surface RAS coordinates of electrodes after 
%           correction for postimplant brain shift. Depth electrode coordinates 
%           are the same as in *.CT 
%     *.PIALVOX: The pial surface voxel coordinates of electrodes after 
%           correction for postimplant brain shift. Voxel coordinates are 
%           for brainmask.nii.gz file also in the elec_recon folder. Depth 
%           electrode coordinates are the same as in *.CT
%     *.electrodeNames: A text file that indicates the name, type of 
%           electrode (strip, grid, depth), and hemisphere in which each 
%           electrode lies. The ith row of this file corresponds to the ith 
%           row of coordinate files.
%    localization_process_date.log - Record of command line output produced
%                                    when this function is run
%
% In the above, *=Freesurfer subject name and
% date=the date on which those files were generated
%
% Note, depth electrode coordinates are not affected by this function. They
% are kept the same as in the postimplant CT or MRI scan.
%
%
% Authors:
% Andrew Dykstra & David M. Groppe
% June 2015

% To Do
% -Create text files of avg brain coordinates too?

if nargin<2,
   minimizeChange=1; 
end

fsDir=getFsurfSubDir();

%% Get version of code to store with electrode coordinates:
brainShiftMethod=['dykstra-' iELVis_getGitInfo];

%%
[elecMatrix, elecLabels, elecRgb, elecPairs, elecPresent]=mgrid2matlab(sub);

% Remove electrodes that were disabled in bioimagesuite
presentIds=find(elecPresent);
elecMatrix=elecMatrix(presentIds,:);
elecMatrix=elecMatrix-1;
elecLabels=elecLabels(presentIds);
elecRgb=elecRgb(presentIds,:);
nElec=length(elecLabels);

elecHem=zeros(nElec,1); %Which hemisphere each electrode is in/over, Left=1, Right=0;
elecStems=cell(nElec,1);
elecNums=zeros(nElec,1);
elecType=cell(nElec,1);
for a=1:nElec,
    if strcmpi(elecLabels{a}(1),'L')
        elecHem(a)=1;
    end
    if ~sum(upper(elecLabels{a}(2))=='GSD')
        error('2nd character of electrode label %s needs to be a G, S, or D.', ...
            elecLabels{a});
    else
        elecType{a}=upper(elecLabels{a}(2));
    end
    underIds=find(elecLabels{a}=='_');
    elecStems{a}=elecLabels{a}(underIds(1)+1:underIds(2)-1);
    elecNums(a)=str2num(elecLabels{a}(underIds(2)+1:end));
end


%% Start diary
elecReconPath=fullfile(fsDir,sub,'elec_recon');

diary_file = fullfile(elecReconPath,['localization_process_' datestr(now,29) '.log']);
fprintf('Recording command line output in file: \n%s\n',diary_file);
diary(diary_file)

fprintf('\n================================================================\n');
fprintf('Starting Dykstra et al. localization process for %s at %s\n',sub,datestr(now,31));
fprintf('Freesurfer Recon dir: %s\n',elecReconPath);
fprintf('Initial location mgrid file: %s.mgrid\n',sub);

% Dykstra code uses RAS coordinates (not VOX)
VOX2RAS=[-1 0 0 128; 0 0 -1 128; 0 -1 0 128; 0 0 0 1];
ctRAS=(VOX2RAS*[elecMatrix'; ones(1, size(elecMatrix,1))])';
ctRAS=ctRAS(:,1:3);

%% Find depth electrodes
isDepth=ones(nElec,1);
for a=1:nElec,
    if lower(elecLabels{a}(2))~='d'
        isDepth(a)=0;
    end
end
depthIds=find(isDepth);
sduralIds=find(~isDepth);
fprintf('%d Depth electrodes\n',length(depthIds));
fprintf('%d Subdural electrodes\n',length(sduralIds));


%% Initialize Coordinate Variables
leptoRAS=ctRAS;
leptoVOX=zeros(nElec,3);
infRAS=zeros(nElec,3)*NaN;
pialRAS=ctRAS;
pialVOX=zeros(nElec,3);


%% Correct subdurals for brain shift
% The optimset function used by Dykstra's code can't deal with initial
% coordinates that have a value of 0. If you add a tiny amount it
% runs. -DG
jit=zeros(size(ctRAS,1),1);
jitCtRAS=ctRAS;
for a=1:3,
    ids=find(ctRAS(:,a)==0);
    if ~isempty(ids),
        jitCtRAS(ids,a)=.000000000000001;
        jit(ids)=1;
    end
end
if sum(jit)
    fprintf('A small value was added to %d electrodes to prevent their CT coordinates from having a zero.\n', ...
        sum(jit));
end

%surfPath=sprintf('/Applications/freesurfer/subjects/%s/surf/',sub);
surfPath=fullfile(fsDir,sub,'surf');
for hemLoop=0:1,
    hemIds=find(elecHem==hemLoop);
    
    if ~isempty(hemIds)
        if hemLoop==0,
            hem='r';
        else
            hem='l';
        end
        
        %% Load Leptomeningeal Surface
        surftype='pial-outer-smoothed';
        [surf.vert surf.tri]=read_surf(fullfile(surfPath,[hem 'h.' surftype]));
        
        %% Brain Shift Correction
        useIds=intersect(hemIds,sduralIds);
        if ~isempty(useIds)
            if universalYes(minimizeChange)
                % Project subdural electrodes to leptomeningeal surface while minimizing
                % distortion
                coordLepto=snap2dural_energy(jitCtRAS(useIds,:),surf);
            else
                %Simply assign subdural electrodes to the nearest leptomeningeal vertex
                coordLepto=get_loc_snap_mgh(jitCtRAS(useIds,:),surfPath,hem,'pial-outer-smoothed');
            end
            leptoRAS(useIds,:)=coordLepto;
            
            %% Project leptomeningeal locations to pial surface
            pialRAS(useIds,:)=get_loc_snap_mgh(coordLepto,surfPath,hem,'pial');
        end
    end
end


%% Save the electrodes locations and labels as text files

%%%%%% Output Electrode Names to Text Files %%%%%%%%%
fnameLabels=fullfile(elecReconPath,[sub '.electrodeNames']);
fprintf('Saving electrode labels to: %s\n',fnameLabels);
fidLabels=fopen(fnameLabels,'w');
fprintf(fidLabels,'%s\n',datestr(now));
fprintf(fidLabels,'Name, Depth/Strip/Grid, Hem\n');
for a=1:nElec,
    if elecHem(a)
        hem='L';
    else
        hem='R';
    end
    fprintf(fidLabels,'%s%d %s %s\n',elecStems{a},elecNums(a),elecType{a},upper(hem(1)));
end
fclose(fidLabels);

%%%%%% Output RAS Coordinates to Text Files %%%%%%%%%
% POSTIMPLANT (CT or MRI) RAS COORDINATES
fnamePostImpRAS = fullfile(elecReconPath,[ sub '.POSTIMPLANT']);
fprintf('Saving CT RAS electrode locations to: %s\n',fnamePostImpRAS);
fidPostImp=writeElecCoordHeader(fnamePostImpRAS,brainShiftMethod,sub);
for a=1:nElec,
    fprintf(fidPostImp,'%f %f %f\n',ctRAS(a,1),ctRAS(a,2),ctRAS(a,3));
end
fclose(fidPostImp);

% Lepto RAS COORDINATES
fnameLeptoRAS = fullfile(elecReconPath,[sub '.LEPTO']);
fprintf('Saving Lepto RAS electrode locations to: %s\n',fnameLeptoRAS);
fidLepto=writeElecCoordHeader(fnameLeptoRAS,brainShiftMethod,sub);
for a=1:nElec,
    fprintf(fidLepto,'%f %f %f\n',leptoRAS(a,1),leptoRAS(a,2),leptoRAS(a,3));
end
fclose(fidLepto);

% Pial RAS COORDINATES
fnamePialRAS = fullfile(elecReconPath,[sub '.PIAL']);
fprintf('Saving Pial RAS electrode locations to: %s\n',fnamePialRAS);
fidPial=writeElecCoordHeader(fnamePialRAS,brainShiftMethod,sub);
for a=1:nElec,
    fprintf(fidPial,'%f %f %f\n',pialRAS(a,1),pialRAS(a,2),pialRAS(a,3));
end
fclose(fidPial);

%%%%%% Output VOX Coordinates to Text Files %%%%%%%%%
% Lepto VOX COORDINATES
RAS2VOX=inv(VOX2RAS);
leptoVOX=(RAS2VOX*[leptoRAS'; ones(1, nElec)])';
fnameLeptoVOX = fullfile(elecReconPath,[sub '.LEPTOVOX']);
fprintf('Saving lepto VOX electrode locations to: %s\n',fnameLeptoVOX);
fidLeptoVox=writeElecCoordHeader(fnameLeptoVOX,brainShiftMethod,sub);
for a=1:nElec,
    fprintf(fidLeptoVox,'%f %f %f\n',leptoVOX(a,1),leptoVOX(a,2),leptoVOX(a,3));
end
fclose(fidLeptoVox);

% Pial VOX COORDINATES
pialVOX=(RAS2VOX*[pialRAS'; ones(1, nElec)])';
fnamePialVOX = fullfile(elecReconPath,[sub '.PIALVOX']);
fprintf('Saving pial VOX electrode locations to: %s\n',fnamePialVOX);
fidPialVox=writeElecCoordHeader(fnamePialVOX,brainShiftMethod,sub);
for a=1:nElec,
    fprintf(fidPialVox,'%f %f %f\n',pialVOX(a,1),pialVOX(a,2),pialVOX(a,3));
end
fclose(fidPialVox);

%% Created text file of Inflated Pial Surface Coordinates (relies on just created text files) 
infRAS=pial2InfBrain(sub,[]);
fnameInfRAS = fullfile(elecReconPath,[sub '.INF']);
fprintf('Saving inflated pial RAS electrode locations to: %s\n',fnameInfRAS);
fidInf=writeElecCoordHeader(fnameInfRAS,brainShiftMethod,sub);
for a=1:nElec,
    fprintf(fidInf,'%f %f %f\n',infRAS(a,1),infRAS(a,2),infRAS(a,3));
end
fclose(fidInf);


%% Plot results to double check
plotPostImpVsLepto(sub,1,1);

% close diary
fprintf('\nElectrodes Localization finished for %s',sub);
fprintf('\n================================================================\n');
diary off


