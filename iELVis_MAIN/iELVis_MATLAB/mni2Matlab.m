function [elecMatrix, elecLabels, elecRgb, elecPairs, elecPresent]=mni2Matlab(sub,elecHem)
%function [elecMatrix, elecLabels, elecRgb, elecPairs, elecPresent]=mni2Matlab(sub,elecHem)
%
% Required Input:
%  sub - Freesurfer name of the subject.
%
% Optional Inputs:
%  elecHem - 'L', 'R', or 'FirstChar': If 'L' all electrodes are assumed to lie on
%        the left hemisphere. If 'R' all electrodes are assumed to lie
%        onthe right hemisphere. 'FirstChar' means that each electrode's
%        name begins with 'L' or 'R', which specifies the hemisphere. If
%        empty, the electrode's hemisphere assignment is automatically
%        deteremined using its anatomical location. Automatic assignment
%        may fail for medial electrodes and can be corrected by manually
%        editing the patient's *.electrodeNames file. Default: elecHem=[]
%
% Output:
%  elecMatrix  - 3D matrix of electrode voxel coordinates in LIP order.
%                Column 1 is R->L (i.e., 1 is the right-most plane of voxels).
%                Column 2 is in S->I (i.e., 1 is the most superior plane of
%                voxels). Column 3 is A->P (i.e., 1 is the most anterior
%                plane of voxels). NOTE THAT 1 IS ADDED TO EACH MGRID
%                COORDINATE, since BioimageSuite starts indexing at 0. This
%                is the same coordinate system output by mgrid2matlab.m.
%  elecLabels  - Cell array of electrode names corresponding to each row of
%                elecMatrix. Each contact's stem is separated from its
%                number by an underscore.
%  elecRgb     - Matrix of RGB colors for each electrode. The nth row of
%                elecRgb corresponds to the nth element of elecLabels.
%  elecPairs   - n x 3 cell array. First two columns indicate electrode
%                pairs that are neighbors. Last column indicates the color
%                (rgb) of the pair. For use with plotElecPial.m.
%  elecPresent - Binary vector of length(elecMatrix). If elecPresent(x)==0
%                then the xth electrode is disabled (i.e., cut off and not
%                implanted).
%
%  Example:
%   >>[elecMatrix, elecLabels, elecRgb]=mni2Matlab('PT001');
%
% Note, this is designed to produce the exact same outputs as
% mgrid2matlab.m. To do this, it has to convert electrode coordinates from
% RAS to LIP. It assumes that the FreeSurfer nii volume dimensions are
% 256^3.
%


fsDir=getFsurfSubDir();
elecReconDir=fullfile(fsDir,sub,'elec_recon');

%setenv('PATH', [getenv('PATH') ':/usr/local/fsl/bin']);
%/usr/local/fsl/bin/flirt
setenv('PATH', [getenv('PATH') ':/usr/local/fsl/bin']);

if nargin<2,
    elecHem=[];
end

if ~isempty(elecHem),
    ids=findStrInCell(elecHem,{'L','R','FirstChar'},0);
    if isempty(ids),
        error('Illegal value of elecHem.');
    end
end

defaultElecType='D';

% Load mni electrode information
mniInfoFname=fullfile(elecReconDir,'mniElecInfo.tsv');
if ~exist(mniInfoFname,'file')
    error('mni file %s not found.',mniInfoFname);
end
infoCsv=csv2Cell(mniInfoFname,9,1);
nElec=size(infoCsv,1);
mniRas=zeros(nElec,3); % RAS coordinates in mni space
elecRgb=zeros(nElec,3); % RGB colors for each electrod
elecLabels=cell(nElec,1); % Electrode names
simpleLabels=cell(nElec,1); % Simplified electrode names like "LACC1" (for use in this function only)
elecPresent=zeros(nElec,1); % 0=electrode cut out
for a=1:nElec,
    % Get Coordinates
    for b=1:3,
        mniRas(a,b)=str2double(infoCsv{a,b+2});
    end
    
    % Elec Present
    if strcmpi(infoCsv{a,6},'TRUE')
        elecPresent(a)=1;
    end
    
    % Simple labels
    simpleLabels{a}=[infoCsv{a,1} infoCsv{a,2}];
    
    %     % Elec Type
    %     if strcmpi(infoCsv{a,11},'Depth'),
    %         elecType='D';
    %     elseif strcmpi(infoCsv{a,10},'TRUE')
    %         elecType='G';
    %     elseif strcmpi(infoCsv{a,11},'Subdural')
    %         elecType='S';
    %     else % might be Unspecified
    %         fprintf('Electrode type for %s is %s. Setting it to %s by default.', ...
    %             elecType,infoCsv{a,1},infoCsv{a,11},defaultElecType);
    %     end
    %     % Elec Name
    %     hem=infoCsv{a,1}(1); % TODO in need to make this a separate column in mniInfo or generate it automatically
    %     elecLabels{a}=sprintf('%c%c_%s_%s',hem,elecType,infoCsv{a,1}, ...
    %         infoCsv{a,2}); % elecLabels are formatted like this LS_LO_5
    % Elec RGB
    for b=1:3,
        elecRgb(a,b)=str2double(infoCsv{a,b+6})/255;
    end
end

%%
if 0,
    % Apply FLIRT transformation within Matlab THIS DOESN'T WORK. FLIRT
    % must get some information from the nii header to change how it
    % applies the matrix.
    % Load matrix mapping from mni to FreeSurfer space
    transFname=fullfile(elecReconDir,'mni2fsurf.mat'); % Affine transformation matrix produced by flirt
    fid = fopen(transFname);
    affineTrans = textscan(fid,'%f%f%f%f','HeaderLines',0,'CollectOutput',1);
    affineTrans = affineTrans{:};
    fclose(fid);
    
    % Convert mni coordinates to Mgrid-style coordinates
    fsurfRas=affineTrans*[mniRas ones(nElec,1)]'; % convert coordinates from 3rd party app to FreeSurfer space
   
else    
    % Import electrode locations derived from FLIRT's img2imgcoord
    temp=importdata(fullfile(elecReconDir,'fsurfXyz.txt'));
    fsurfRas=temp.data'; % coordinates are Right+, Ant+, Sup+ 
end
elecMatrix=zeros(nElec,3); % LIP coordinates
% Note, fsurfRas coordinates range from 0-255, so elecMatrix coordinates
% will range from 1-256 which is exactly the range produced by
% mgrid2matlab.m
% For some reason this works for TWH cases
elecMatrix(:,1)=-1+fsurfRas(1,:); % convert R to L %elecMatrix(:,1)=256-fsurfRas(1,:); % convert R to L
elecMatrix(:,2)=1+fsurfRas(2,:); % it's already I, but origin needs to be 1
elecMatrix(:,3)=256-fsurfRas(3,:); % convert A to P
% But this works for PT001 and makes more sense
% elecMatrix(:,1)=256-fsurfRas(1,:); % convert R to L
% elecMatrix(:,2)=256-fsurfRas(3,:); % convert S to I
% elecMatrix(:,3)=256-fsurfRas(2,:); % convert A to P

%% For debugging
% figure(1); clf;
% subplot(1,2,1);
% plot3(fsurfRas(1,:),fsurfRas(2,:),fsurfRas(3,:),'.','markersize',14);
% xlabel('ONE');
% ylabel('TWO');
% zlabel('THREE');
% title('fsurfRAS');
% axis square;
% 
% subplot(1,2,2);
% plot3(elecMatrix(:,1),elecMatrix(:,2),elecMatrix(:,3),'.','markersize',14);
% xlabel('ONE-L');
% ylabel('TWO-I');
% zlabel('THREE-P');
% title('elecMatrix');
% axis square;

%% Get elec labels
elecCoordILA=zeros(nElec,3);
elecCoordILA(:,1)=elecMatrix(:,2);
elecCoordILA(:,2)=elecMatrix(:,1);
elecCoordILA(:,3)=fsurfRas(2,:);
for a=1:nElec,
    fprintf('Identifying hemisphere corresponding to electrode %d/%d\n',a,nElec);
    hem=[];
    if strcmpi(elecHem,'L'),
        hem='L';
    elseif strcmpi(elecHem,'R'),
        hem='R';
    elseif strcmpi(elecHem,'FirstChar'),
        hem=infoCsv{a,1}(1);
    else
        anatLabel=vox2Seg(round(elecCoordILA(a,:)),'PT001');
        seg=strsplit(anatLabel,'-');
        
        if length(seg)>1,
            if strcmpi(seg(2),'lh'),
                hem='L';
            elseif strcmpi(seg(2),'rh'),
                hem='R';
            end
        end
        if isempty(hem),
            % anatomical location unknown, guess left hemisphere if coordinate
            % is on left side of volume (volume dimensions are 256^3)
            if  elecCoordILA(a,2)<=128,
                hem='R';
            else
                hem='L';
            end
        end
    end
    
    % Elec Type
    if strcmpi(infoCsv{a,11},'Depth'),
        elecType='D'; % depth
    elseif strcmpi(infoCsv{a,10},'TRUE')
        elecType='G'; % grid
    elseif strcmpi(infoCsv{a,11},'Subdural')
        elecType='S'; % strip
    else % might be Unspecified
        fprintf('Electrode type for %s is %s. Setting it to %s by default.', ...
            elecType,infoCsv{a,1},infoCsv{a,11},defaultElecType);
    end
    % Elec Name
    elecLabels{a}=sprintf('%c%c_%s_%s',hem,elecType,infoCsv{a,1}, ...
        infoCsv{a,2}); % elecLabels are formatted like this LS_LO_5
end

%%
% Load mni electrode pair information
mniPairsFname=fullfile(elecReconDir,'mniElecPairs.tsv');
if ~exist(mniPairsFname,'file')
    error('mni file %s not found.',mniPairsFname);
end
pairsCsv=csv2Cell(mniPairsFname,9,1);
nPairs=size(pairsCsv,1);
elecPairs=cell(nPairs,3);
for a=1:nPairs,
    % Find the index and iELVis name of the first electrode
    tempLabel=[pairsCsv{a,1} pairsCsv{a,5}];
    eId1=findStrInCell(tempLabel,simpleLabels,1);
    elecPairs{a,1}=elecLabels{eId1};
    
    % Find the index and iELVis name of the second electrode
    tempLabel=[pairsCsv{a,1} pairsCsv{a,6}];
    eId2=findStrInCell(tempLabel,simpleLabels,1);
    elecPairs{a,2}=elecLabels{eId2};
    
    tempRgb=zeros(1,3);
    for b=1:3,
        tempRgb(b)=str2double(pairsCsv{a,b+1})/255;
    end
    elecPairs{a,3}=tempRgb;
end

%TODO scan mgrid2matlab to make sure I'm not missing anything here that I
%should be doing
