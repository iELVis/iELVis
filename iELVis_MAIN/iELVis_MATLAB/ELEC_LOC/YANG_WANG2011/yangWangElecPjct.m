function yangWangElecPjct(sub)
%function yangWangElecPjct(sub)
%
% Corrects intracranial electrode locations for brain shift using the
% following method:
%  Yang, A. I., Wang, X., Doyle, W. K., Halgren, E., Carlson, C., Belcher, 
%  T. L., et al. (2012). Localization of dense intracranial electrode arrays 
%  using magnetic resonance imaging. NeuroImage 63(1), 157?165. 
%  doi:10.1016/j.neuroimage.2012.06.039
%
% Input:
%  sub - Freesurfer subject name (e.g., 'TWH001')
%
% Outputs:
%  The following files are created in the elec_recon subfolder of the
%  Freesufer subject folder:
%    *.PIAL - RAS coordinates snapped to pial surface
%    *.PIALVOX - Voxel coordinates snapped to pial surface
%    *.LEPTO - RAS coordinates snapped to leptomeningeal (i.e., smoothed pial)
%                  surface)
%    *.electrodeNames - electrode names
%    localization_process_date.log - Record of command line output produced
%                                    when this function is run
%
% In the above, *=Freesurfer subject name and
% date=the date on which those files were generated
%
% Note, depth electrode coordinates are not affected by this function. They
% are kept the same as in the postimplant CT or MRI scan.
%
% Also note that you should double check grid labels as there may be a
% mismatch between NYU interpolation and mgrid file.
%
% Authors:
% Hugh Wang (NYU) & David M. Groppe (University of Toronto)
% (Bugs are DG's fault)
% June 2015


% Future work:
% -get mapping to MNI brain to work
% -Hugh's original code had an option to import electrode info from an excel
% file, but I have not looked to see if it works.
% -Hugh's original code could interpolate depth electrode locations. I've
% disabled it but you could make it an option -DG

% get the subject info
fsDir=getFsurfSubDir();

subPath = fullfile(fsDir,sub);
elecReconPath=fullfile(subPath,'elec_recon');

if ~isdir(subPath)
    error('Freesurfer folder %s not found',subPath);
end


%% read the inital text file with postimplant electrode coordinates
% Note coordinates are in voxels (NOT RAS)
postimpLocFname=fullfile(elecReconPath,[sub 'PostimpLoc.txt']);
if ~exist(postimpLocFname,'file')
    error('File %s does NOT exist.',postimpLocFname);
end


%% Map to MNI brain
% The code below is necessary for doing mapping to MNI brain, DG commented
% them out for the time being
% [dep_img_file,dep_img_path] = uigetfile({'*.nii.gz';'*.nii';'*.img'},'Select the T1 pre-operation image: ',PathName);
% dep_img_path='/Applications/freesurfer/subjects/TWH008/elec_recon/';
% dep_img_file='T1.nii.gz';

%% Start diary
diary_file=fullfile(elecReconPath,['localization_process_' datestr(now,29) '.log']);
fprintf('Recording command line output in file: \n%s\n',diary_file);
diary(diary_file)

fprintf('\n================================================================\n');
fprintf('Starting Yang, Wang, et al. localization process for %s at %s\n',sub,datestr(now,31));
fprintf('Freesurfer Recon dir: %s\n',elecReconPath);
fprintf('Initial location text file: %s\n',postimpLocFname);


%% Read in preop T1
t1Fname=fullfile(elecReconPath,'T1.nii.gz');
hdr = ntools_elec_load_nifti(t1Fname,1);
if ~isequal(hdr.pixdim(2),hdr.pixdim(3),hdr.pixdim(4))
    warning('T1 voxel mm dimensions not equal. This will affect the accuracy of distance calculation');
end
scale = mean(hdr.pixdim(2:4));


%% Read in postimp elec locations text/xls file
[~,~,ext] = fileparts(postimpLocFname);

if strcmpi(ext,'.txt')
    fid = fopen(postimpLocFname);
    elec_all = textscan(fid,'%s %d %f %f %f %s %s','CommentStyle','%');
    elec_cell = [elec_all{1},num2cell(elec_all{2}),num2cell(elec_all{3}), ...
        num2cell(elec_all{4}),num2cell(elec_all{5}),elec_all{6},elec_all{7}];
    fclose(fid);
    % Note the order of the info in elec_cell is:
    % 1) Chan stem
    % 2) Chan #
    % 3-5) Chan voxel coordinates
    % 6) Chan hemisphere
    % 7) Electrode type (D, S, or G)
elseif strcmpi(ext,'.xls') || strcmpi(ext,'.xlsx')
    % I haven't tested this to make sure that it works. DG
    [~,~,elec_cell] = xlsread(postimpLocFname);
end

% split elecs by type
g = strmatch('G',upper(elec_cell(:,7)));
d = strmatch('D',upper(elec_cell(:,7)));
s = strmatch('S',upper(elec_cell(:,7)));
L = strmatch('L',upper(elec_cell(:,6)));
R = strmatch('R',upper(elec_cell(:,6)));

nManGrid=length(g);
nStrip=length(s);
nDepth=length(d);
nManElec=nDepth+nManGrid+nStrip;
nHem=[0 0];
nHem(1)=length(L);
nHem(2)=length(R);

fprintf('%d Depth electrodes manually marked\n',nDepth);
fprintf('%d Grid electrodes manually marked\n',nManGrid);
fprintf('%d Strip electrodes manually marked\n',nStrip);
fprintf('%d left hemisphere electrodes \n',nHem(1));
fprintf('%d right hemisphere electrodes\n',nHem(2));

% Convert all coordinates to RAS
postimpVox=zeros(nManElec,3);
fullElecNames=cell(nManElec,1);
for a=1:nManElec,
    for b=1:3,
        postimpVox(a,b)=elec_cell{a,b+2};
    end
    fullElecNames{a}=[elec_cell{a,1} num2str(elec_cell{a,2})];
end
VOX2RAS=[-1 0 0 128; 0 0 -1 128; 0 -1 0 128; 0 0 0 1];
postimpRas=(VOX2RAS*[postimpVox'; ones(1, nManElec)])';
%fprintf('RAS coordinates:\n');
for a=1:nManElec,
    for b=1:3,
        elec_cell{a,b+2}=postimpRas(a,b);
    end
    %fprintf('%s: %.2f %.2f %.2f\n',elec_cell{a,1},elec_cell{a,2}, ...
    %    elec_cell{a,3},elec_cell{a,4});
end

% outer-brain surface check and create
ntools_elec_outer_brain(subPath); % If smoothed pial surface has NOT been created, I believe this function fails. DG

% Initialize text file Ids
fidLepto=[];
fidPial=[];
fidCT=[];
fidLabels=[];
fidLeptoVox=[];
fidPialVox=[];

%% Calculate electrode locations with correction for brain shift
%% Loop over hemispheres
for hemLoop=1:2,
    if hemLoop==1
        hem='lh';
        hemIds=L;
    else
        hem='rh';
        hemIds=R;
    end
    
    if nHem(hemLoop),
        elec_depth=elec_cell(intersect(d,hemIds),:);
        elec_depthCT=elec_depth;
        elec_gridCT=elec_cell(intersect(g,hemIds),:);
        elec_strip=elec_cell(intersect(s,hemIds),:);
        elec_stripCT=elec_strip;
        %fullElecNamesThisHem=fullElecNames(hemIds);
        nGridThisHem=size(elec_gridCT,1);
        nDepthThisHem=size(elec_depthCT,1);
        nStripThisHem=size(elec_stripCT,1);
               
        if nGridThisHem,
            % Figure out how many unique types of grids there are
            gridNames = regexp(elec_gridCT(:,1),'[A-Za-z]*[^\d*]','match');
            ini_gridNames=cell(length(gridNames),1);
            for i=1:length(gridNames)
                ini_gridNames(i) = gridNames{i};
            end
            gridNames = unique(ini_gridNames); %# of unique grids
            nGridType=length(gridNames);
            elec_grid=[];
            grid_stats=[];
            
            for a=1:nGridType,
                % For each grid identify the dimensions
                disp('If grids have unequal dimensions, rows is probably the smaller dimension.');
                nRow=input(sprintf('How many rows does %s have? (default-> 8): ',gridNames{a}));
                if isempty(nRow),
                    nRow=8;
                end
                nCol=input(sprintf('How many columns does %s have? (default-> 8): ',gridNames{a}));
                if isempty(nCol),
                    nCol=8;
                end
                
                % For each grid identify the corners in clockwise or counter clockwise
                % order
                defaultCorners=[1 nCol nCol*nRow nCol*nRow-nCol+1];
                corners=input(sprintf('Enter the electrode #''s of %s''s corners starting at 1 and going counterclockwise with square brackets (default-> [%d %d %d %d])', ...
                    gridNames{a},defaultCorners(1),defaultCorners(2),defaultCorners(3),defaultCorners(4)));
                if isempty(corners),
                    %corners=[1 8 64 57];
                    corners=defaultCorners;
                end
                % Select out the corner electrodes to pass to Hugh's function:
                cornerIds=zeros(4,1);
                for b=1:4,
                    %cornerIds(b)=findStrInCell([gridNames{a} int2str(corners(b))],elec_gridCT(:,1),1);
                    cornerIds(b)=findStrInCell([gridNames{a} int2str(corners(b))],fullElecNames(intersect(g,hemIds)),1);
                end
                
                %radius = input('\nInput the inter-electrode distance (mm) (Press Enter for default distance 10mm) : ');
                %if isempty(radius)
                %    radius = 10;
                %end
                radius=10;
                [elec_gridTemp, grid_statsTemp] = ntoolsElecCalcGrid(elec_gridCT(cornerIds,:),subPath,scale,radius,nRow,nCol,hem);
                
                % Collect all temp grid coordinates
                elec_grid=[elec_grid; elec_gridTemp];
            end
        end
        %nGrid=size(elec_grid,1); % Right now, all grid elecs must be labeled in CT
        %scan. If we change this so that only corners are labelled. We have to
        %recount the number of grid elecs here. DG
        nElecThisHem=nDepthThisHem+nGridThisHem+nStripThisHem;
        
        %% Process depth elecs
        %         if 1
        %             fprintf('Localizing depths using manually marked locations.\n');
        %             fprintf('Locations are NOT corrected for brain shift.\n');
        %         else
        %             % Hugh's original code
        %             fprintf('Localizing depths using the most extreme electrodes and interpolating the rest.\n');
        %             elec_depth = ntools_elec_calc_depth(ini_depth);
        %         end
        
        %% Process strip elecs
        if nStripThisHem,
            elec_strip = ntools_elec_calc_strip(elec_strip,subPath,hem);
        end
        
        %% Snap subdural electrodes to the pial surface
        ct=0;
        nSnapThisHem=nGridThisHem+nStripThisHem;
        pialRAS=zeros(nSnapThisHem+nDepthThisHem,3);
        for a=1:nGridThisHem,
            ct=ct+1;
            for b=1:3,
                pialRAS(ct,b)=elec_grid{a,b+2};
            end
        end
        for a=1:nStripThisHem,
            ct=ct+1;
            for b=1:3,
                pialRAS(ct,b)=elec_strip{a,b+2};
            end
        end
        if nSnapThisHem,
            pialRAS = snap2surf(pialRAS(1:ct,:),fullfile(subPath,'surf'),hem(1),'pial');
        end
        
        %% Add depth coords, which are not snapped to surface
        for a=1:nDepthThisHem,
            ct=ct+1;
            for b=1:3,
                pialRAS(ct,b)=elec_depth{a,b+2};
            end
        end
        
        
        %% Save the electrodes locations as text files(note that the text 
        % files are closed after looping over both hemispheres)
        leptoRAS=zeros(nElecThisHem,3);
        ctRAS=leptoRAS;
        elecStems=cell(nElecThisHem,1);
        elecNums=zeros(nElecThisHem,1);
        elecType=cell(nElecThisHem,1);
        ct=0;
        for a=1:nGridThisHem, % ?? make sure this works
            ct=ct+1;
            elecStems{ct}=elec_grid{a,1};
            elecNums(ct)=elec_grid{a,2};
            elecType{ct}='G';
            tempIdStem=findStrInCell(elecStems{ct},elec_gridCT(:,1),1);
            tempIdNum=find(cellfun(@(x) isequal(x,elecNums(ct)),elec_gridCT(:,2)));
            tempId=intersect(tempIdNum,tempIdStem);
            if isempty(tempId)
               error('Could not find channel %s%d in elec_gridCT',elecStems{ct}, ...
                   elecNums(ct));
            end
            for b=1:3,
                leptoRAS(ct,b)=elec_grid{a,b+2};
                ctRAS(ct,b)=elec_gridCT{tempId,b+2};
            end
        end
        for a=1:nStripThisHem,
            ct=ct+1;
            for b=1:3,
                leptoRAS(ct,b)=elec_strip{a,b+2};
                ctRAS(ct,b)=elec_stripCT{a,b+2};
            end
            elecStems{ct}=elec_strip{a,1};
            elecNums(ct)=elec_strip{a,2};
            elecType{ct}='S';
        end
        for a=1:nDepthThisHem,
            ct=ct+1;
            for b=1:3,
                leptoRAS(ct,b)=elec_depth{a,b+2};
                ctRAS(ct,b)=elec_depthCT{a,b+2};
            end
            elecStems{ct}=elec_depth{a,1};
            elecNums(ct)=elec_depth{a,2};
            elecType{ct}='D';
        end
        
        % RAS COORDINATES
        % Lepto
        fnameLeptoRAS = fullfile(elecReconPath,[sub '.LEPTO']);
        fprintf('Saving lepto RAS electrode locations to: %s\n',fnameLeptoRAS);
        if isempty(fidLepto)
            fidLepto=fopen(fnameLeptoRAS,'w');
            fprintf(fidLepto,'%s\n',datestr(now));
            fprintf(fidLepto,'R A S\n');
        end
        for a=1:nElecThisHem,
            fprintf(fidLepto,'%f %f %f\n',leptoRAS(a,1),leptoRAS(a,2),leptoRAS(a,3));
        end
        
        % Pial
        fnamePialRAS = fullfile(elecReconPath,[sub '.PIAL']);
        fprintf('Saving pial RAS electrode locations to: %s\n',fnamePialRAS);
        if isempty(fidPial)
            fidPial=fopen(fnamePialRAS,'w');
            fprintf(fidPial,'%s\n',datestr(now));
            fprintf(fidPial,'R A S\n');
        end
        for a=1:nElecThisHem,
            fprintf(fidPial,'%f %f %f\n',pialRAS(a,1),pialRAS(a,2),pialRAS(a,3));
        end
        
        % CT (i.e., uncorrected for brain shift)
        fnameCtRAS = fullfile(elecReconPath,[sub '.CT']);
        fprintf('Saving CT RAS electrode locations to: %s\n',fnameCtRAS);
        if isempty(fidCT)
            fidCT=fopen(fnameCtRAS,'w');
            fprintf(fidCT,'%s\n',datestr(now));
            fprintf(fidCT,'R A S\n');
        end
        for a=1:nElecThisHem,
            fprintf(fidCT,'%f %f %f\n',ctRAS(a,1),ctRAS(a,2),ctRAS(a,3));
        end
        
        % Electrode names 
        fnameLabels = fullfile(elecReconPath,[sub '.electrodeNames']);
        fprintf('Saving electrode labels to: %s\n',fnameLabels);
        if isempty(fidLabels)
            fidLabels=fopen(fnameLabels,'w');
            fprintf(fidLabels,'%s\n',datestr(now));
            fprintf(fidLabels,'Name, Depth/Strip/Grid, Hem\n');
        end
        for a=1:nElecThisHem,
            fprintf(fidLabels,'%s%d %s %s\n',elecStems{a},elecNums(a),elecType{a},upper(hem(1)));
        end
        
        % VOX COORDINATES
        RAS2VOX=inv(VOX2RAS);
        leptoVOX=(RAS2VOX*[leptoRAS'; ones(1, nElecThisHem)])';
        fnameLeptoVOX = fullfile(elecReconPath,[sub '.LEPTOVOX']);
        fprintf('Saving lepto VOX electrode locations to: %s\n',fnameLeptoVOX);
        if isempty(fidLeptoVox)
            fidLeptoVox=fopen(fnameLeptoVOX,'w');
            fprintf(fidLeptoVox,'%s\n',datestr(now));
            fprintf(fidLeptoVox,'X Y Z\n');
        end
        for a=1:nElecThisHem,
            fprintf(fidLeptoVox,'%f %f %f\n',leptoVOX(a,1),leptoVOX(a,2),leptoVOX(a,3));
        end
        
        pialVOX=(RAS2VOX*[pialRAS'; ones(1, nElecThisHem)])';
        fnamePialVOX = fullfile(elecReconPath,[sub '.PIALVOX']);
        fprintf('Saving pial VOX electrode locations to: %s\n',fnamePialVOX);
        if isempty(fidPialVox)
            fidPialVox=fopen(fnamePialVOX,'w');
            fprintf(fidPialVox,'%s\n',datestr(now));
            fprintf(fidPialVox,'X Y Z\n');
        end
        for a=1:nElecThisHem,
            fprintf(fidPialVox,'%f %f %f\n',pialVOX(a,1),pialVOX(a,2),pialVOX(a,3));
        end
        
    end
end

%% Close text files
fclose(fidLepto);
fclose(fidPial);
fclose(fidCT);
fclose(fidLabels);
fclose(fidLeptoVox);
fclose(fidPialVox);


%% Created text file of Inflated Pial Surface Coordinates (relies on recently closed text files) 
infRAS=pial2InfBrain(sub,[]);
fnameInfRAS = fullfile(elecReconPath,[sub '.INF']);
fprintf('Saving inflated pial RAS electrode locations to: %s\n',fnameInfRAS);
fidInf=fopen(fnameInfRAS,'w');
fprintf(fidInf,'%s\n',datestr(now));
fprintf(fidInf,'R A S\n');
for a=1:nManElec,
    fprintf(fidInf,'%f %f %f\n',infRAS(a,1),infRAS(a,2),infRAS(a,3));
end
fclose(fidInf);


%% Plot results to double check
plotCtVsLepto(sub,1,1);


%% save all into binary nifti image
if 0
    % DG disabled this as I don't know how useful it is
    %fname_bin = [PathName,Sname,'_elec_bin_T1_' datestr(now,29), '.nii.gz'];
    fname_bin = fullfile(elecReconPath,[sub '_elec_bin_T1_' datestr(now,29), '.nii.gz']);
    fprintf('Saving electrode locations as binary nii file: %s\n',fname_bin);
    elec_vox = ntools_elec_savebin([x y z],hdr,fname_bin);
end

% transform into MNI space  ?? make work later
if 0
    elec_mni = ntools_elec_dartel_warp(fname_bin,[dep_img_path,dep_img_file]);
    %fname_mni = [PathName Sname '_coor_MNI_' datestr(now,29) '.txt'];
    fname_mni = fullfile(elecReconPath,[sub '_coor_MNI_' datestr(now,29) '.txt']);
    ntools_elec_savetxt(fname_mni,[name num2cell(elec_mni) label]);
end


% close diary
fprintf('\nElectrodes Localization finished for %s',sub);
fprintf('\n================================================================\n');
diary off


end

