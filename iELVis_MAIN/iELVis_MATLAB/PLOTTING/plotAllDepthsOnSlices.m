function plotAllDepthsOnSlices(fsSub,elecInfoType,cfg)
% function plotAllDepthsOnSlices(fsSub,elecInfoType,cfg)
%
% Creates a figure illustrating the location of each depth electrode contact 
% in a sagittal, coronal, and axial slice and indicates which part of
% the brain it is in. You need to first need to record electrode
% coordinates in iELVis conventions (e.g., using yangWangElecPjct.m or
% dykstraElecPjct.m), have an mgrid file of electrode information from
% BioImageSuite, or have analogous info in MNI space.
%
% Required Inputs:
%  fsSub - Patient's freesurfer directory name
%  elecInfoType - {'mgrid', 'BIDS-iEEG', or 'MNI'} a string indicating the
%                 format electrode information is stored in. 
%
% Optional cfg parameters:
%  mgridFname - mgrid filename and path. If empty, and
%               elecInfoType=='mgrid', we assume the mgrid file is in the
%               subject's elec_recon subfolder and named *.mgrid, where *
%               is the subject's FreeSurfer ID.
%  fullTitle  - If 1, the mgrid and mri voxel coordinates are displayed in
%               the figure title along with the electrode name and anatomical
%               location. {default: 0}
%  markerSize - The size of the dot in each slice used to represent an
%               electrode's location. {default: 30}
%  cntrst     - 0< number <=1 The lower this number, the lower the brightness
%               of the image (i.e., the lower the voxel value corresponding to
%               white). {default: 0.5}
%  anatOverlay -If 1, color is overlayed on the brain to show FreeSurfer's
%              automatic segmentation of brain areas (neocortex uses
%              Desikan-Killiany parcellation). Alternatively define the fullpath
%              to another parcellation file. {default: 0}
%  colorLUT    - fullpath to color lookup table if you would like to use
%                non default colors for your parcellation overlay.
%                {default: FreeSurferColorLUTnoFormat.txt}
%  pauseOn   - If 1, Matlab pauses after each figure is made and waits for
%              a keypress. {default: 0}
%  figOverwrite - If non-zero, only one new Matlab figure will be produced 
%              when this function is called. Each time a new electrode
%              needs to be visualized the figure is cleared. This is useful
%              when lots of depths have been used and you're printing the
%              figures.
%  printFigs - 1 or directory. If 1, each figure is output to a jpg file in the patient's
%              elec_recon/PICS folder and the figure is closed after the
%              jpg is created. This is particularly useful for implants with
%              a large number of depth contacts. If a directory, the jpg files
%              get saved there instead. {default: 0}
%  bidsDir   - Full path to an BIDS-iEEG root directory. If specified,
%              electrode location and pial surface files will be
%              imported from this directory. If not specified,
%              these data will be imported from the subject's FreeSurfer directory.
%              Note that all electrodes are colored red when you use this
%              option as BIDS-iEEG does not associate a default color with
%              each electrode.
%  bidsSes   - integer. The BIDS-iEEG session number. This has no effect if
%              bidsDir not specified {default: 1}
%
%
% Examples:
%  %Specify mgrid file and do NOT print
%  cfg=[];
%  cfg.mgridFname='/Applications/freesurfer/subjects/PT001/elec_recon/PT001.mgrid';
%  plotAllDepthsOnSlices('PT001','mgrid',cfg);
%
%  %Use FreeSurfer file structure and print
%  cfg=[];
%  cfg.printFigs=1;
%  plotAllDepthsOnSlices('PT001','mgrid',cfg);
%
%
% Author: David M. Groppe
% Feb. 2015
% Feinstein Institute for Medical Research/Univ. of Toronto


if ~isfield(cfg,'mgridFname'),    mgridFname=[];    else mgridFname=cfg.mgridFname; end
if ~isfield(cfg,'fullTitle'),     fullTitle=0;      else fullTitle=cfg.fullTitle; end
if ~isfield(cfg,'markerSize'),    markerSize=30;    else markerSize=cfg.markerSize; end
if ~isfield(cfg,'cntrst'),    cntrst=.5;          else cntrst=cfg.cntrst; end
if ~isfield(cfg,'anatOverlay'),    anatOverlay=.5;          else anatOverlay=1; end
if ~isfield(cfg,'colorLUT'),    colorLUT=0;          else colorLUT=cfg.colorLUT; end
if ~isfield(cfg,'pauseOn'),    pauseOn=0;          else pauseOn=cfg.pauseOn; end
if ~isfield(cfg,'printFigs'),    printFigs=0;          else printFigs=cfg.printFigs; end
if ~isfield(cfg, 'bidsDir'),      bidsDir=[];         else bidsDir=cfg.bidsDir; end
if ~isfield(cfg, 'bidsSes'),      bidsSes=1;         else bidsSes=cfg.bidsSes; end
if ~isfield(cfg, 'figOverwrite'),  figOverwrite=1;  else figOverwrite=cfg.figOverwrite; end
checkCfg(cfg,'plotAllDepthsOnSlices.m');

% Check for valid arguments
if (nargin<3),
   error('plotAllDepthsOnSlices requires 2 arguments'); 
end
validFormats={'mgrid', 'BIDS-iEEG', 'MNI'};
if isempty(findStrInCell(elecInfoType,validFormats))
    error('%s is not a valid value for elecInfoType parameter',elecInfoType);
end
if strcmpi(elecInfoType,'BIDS-iEEG'),
   if isempty(bidsDir),
      error('You need to specifiy cfg.bidsDir when reading electrode information in BIDS-iEEG format.'); 
   end
end

% FreeSurfer Subject Directory
fsdir=getFsurfSubDir();

% Define path to neuroimaging files
if isempty(bidsDir)
    % Look in FreeSurfer directories
    mriDir=fullfile(fsdir,fsSub,'mri');
else
    % Look in BIDS-iEEG directory
    mriDir=fullfile(bidsDir,'derivatives','iELVis',['sub-' fsSub],'mri');
end

% Load MRI
mriFname=fullfile(mriDir,'brainmask.mgz');
if ~exist(mriFname,'file')
    error('File %s not found.',mriFname);
end
mri=MRIread(mriFname);
%mri.vol is ILA (i.e., S->I, R->L, P->A)
mx=max(max(max(mri.vol)))*cntrst;
mn=min(min(min(mri.vol)));
sVol=size(mri.vol);

% Get electrode info
if strcmpi(elecInfoType,'BIDS-iEEG'),
    % BIDS-iEEG format
    % Get electrode information from iEEG-BIDs directory and transform
    % coordinates to be compatible with those in mgrid files
    %/Users/davidgroppe/Desktop/HandMotor/sub-PT001/ieeg/sub-PT001_ses-01_space-inf_electrodes.tsv
    coordFname=fullfile(bidsDir,['sub-' fsSub],'ieeg',sprintf('sub-%s_ses-%.2d_space-postimplant_electrodes.tsv',fsSub,bidsSes));
    fprintf('Taking electrode info from %s.\n',coordFname);
    
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
    % convert RAS coordinates to LIP
    elecMatrix=zeros(nElecTotal,3);
    elecMatrix(:,1)=-RAS_coor(:,1)+129;
    elecMatrix(:,2)=-RAS_coor(:,3)+129;
    elecMatrix(:,3)=-RAS_coor(:,2)+129;
    
    %% Collect electrode name, type, & hemisphere into elecLabels
    % also define electrode colors (all red for the time being)
    elecLabels=cell(nElecTotal,3); % a cell array of things like LD_LDAm10
    nameId=findStrInCell('name',elecCoordCsv(1,:),1);
    typeId=findStrInCell('type',elecCoordCsv(1,:),1);  % grid, depth, strip
    hemId=findStrInCell('hemisphere',elecCoordCsv(1,:),1); % L,R
    elecRgb=zeros(nElecTotal,3);
    elecRgb(:,1)=1; % make all electrodes red
    for csvLoopA=1:nElecTotal,
        tempName=elecCoordCsv{csvLoopA+1,nameId}; % Name
        switch elecCoordCsv{csvLoopA+1,typeId}
            case 'grid'
                tempType='G';
            case 'strip'
                tempType='S';
            case 'depth'
                tempType='D';
            otherwise
                error('Unrecognized value of electrode "type": %s',elecCoordCsv{csvLoopA+1,typeId});
        end
        tempHem=elecCoordCsv{csvLoopA+1,hemId}; % Hemisphere
        elecLabels{csvLoopA}=sprintf('%s%s_%s',tempHem,tempType,tempName);
    end
elseif strcmpi(elecInfoType,'MNI'),
    % MNI format
    elecReconDir=fullfile(fsdir,fsSub,'elec_recon');
    mniInfoFname=fullfile(elecReconDir,'mniElecInfo.tsv');
    mniPairsFname=fullfile(elecReconDir,'mniElecPairs.tsv');
    if ~exist(mniInfoFname,'file')
        error('Missing %s',mniInfoFname);
    end
    if ~exist(mniPairsFname,'file')
       error('Missing %s', mniPairsFname);
    end
    [elecMatrix, elecLabels, elecRgb]=mni2Matlab(fsSub);
    % elecMatrix coords are LIP
    % elecLabels is a cell array of things like LD_LDAm_10
    % elecRgb is nElec x 3 matrix of 0-1 RGB values
else
    % mgrid format
    if isempty(mgridFname)
        [elecMatrix, elecLabels, elecRgb]=mgrid2matlab(fsSub);
    else
        [elecMatrix, elecLabels, elecRgb]=mgrid2matlab(fsSub,mgridFname);
    end
end
            
nElec=length(elecLabels);
elecMatrix=round(elecMatrix);
xyz=zeros(size(elecMatrix));
xyz(:,1)=elecMatrix(:,2);
xyz(:,2)=elecMatrix(:,1);
xyz(:,3)=sVol(3)-elecMatrix(:,3);

depthElecs=zeros(nElec,1);
for a=1:nElec,
    if strcmpi(elecLabels{a}(2),'D')
        depthElecs(a)=1;
    end
end

if universalYes(anatOverlay)
    % Load segmentation
    if ischar(cfg.anatOverlay)
        segFname = cfg.anatOverlay;
    else
        segFname=fullfile(mriDir,'aparc+aseg.mgz');
        if ~exist(mriFname,'file')
            error('File %s not found.',mriFname);
        end
    end
    seg=MRIread(segFname);
    
    % Load segmentation color table
    if universalNo(colorLUT)
        pathstr = fileparts(which('mgrid2matlab'));
        inFile=fullfile(pathstr,'FreeSurferColorLUTnoFormat.txt');
        if ~exist(inFile,'file')
            error('Could not find file %s',inFile);
        end
    elseif exist(colorLUT,'file')
        inFile = colorLUT;
    else
        error('The defined color lookup table was not found');
    end
    fid=fopen(inFile,'r');
    %fid=fopen('/Applications/freesurfer/FreeSurferColorLUTnoFormat.txt','r');
    tbl=textscan(fid,'%d%s%d%d%d%d');
    fclose(fid);
end

figId=0;
for elecId=1:nElec,
    if depthElecs(elecId)
        if figId<1,
            figId=figure();
        else
            if universalYes(figOverwrite),
                 figure(figId);
            else
                figId=figure();
            end
        end
        clf();
        set(figId,'position',[78 551 960 346],'paperpositionmode','auto');
        
        hm=zeros(1,3);
        figure(figId); clf;
        colormap gray;
        %subplot(131);
        wdth=.35;
        wDelt=.33;
        xStart=-.005;
        yStart=.03;
        ht=.9;
        axes('position',[xStart yStart wdth ht]);
        imagesc(squeeze(mri.vol(:,xyz(elecId,2),:)),[mn mx]);
        axis square;
        set(gca,'xdir','reverse');
        hold on;
        
        if universalYes(anatOverlay)
            % Plot segmentation
            for a=1:sVol(1),
                for b=1:sVol(3),
                    if seg.vol(a,xyz(elecId,2),b)
                        segId=find(tbl{1}==seg.vol(a,xyz(elecId,2),b));
                        tempRgb=double([tbl{3}(segId) tbl{4}(segId) tbl{5}(segId)])/255;
                        hM=patch([-.5 .5 .5 -.5]+b,[-.5 -.5 .5 .5]+a,tempRgb);
                        set(hM,'LineStyle','none','FaceAlpha',0.3);
                    end
                end
            end
        end
        
        % Plot electrode
        hm(1)=plot(xyz(elecId,3),xyz(elecId,1),'r.');
        set(hm(1),'color',elecRgb(elecId,:),'markersize',markerSize);
        %find image limits
        mxX=max(squeeze(mri.vol(:,xyz(elecId,2),:)),[],2);
        mxY=max(squeeze(mri.vol(:,xyz(elecId,2),:)),[],1);
        limXa=max(intersect(1:(sVol(3)/2),find(mxX==0)));
        limXb=min(intersect((sVol(3)/2:sVol(3)),find(mxX==0)));
        limYa=max(intersect(1:(sVol(1)/2),find(mxY==0)));
        limYb=min(intersect((sVol(1)/2:sVol(1)),find(mxY==0)));
        %keep image square
        tempMin=min([limXa limYa]);
        tempMax=max([limXb limYb]);
        if tempMin<tempMax,
            axis([tempMin tempMax tempMin tempMax]);
        end
        set(gca,'xtick',[],'ytick',[]);
        
        %subplot(132);
        axes('position',[xStart+wDelt yStart wdth ht]);
        imagesc(squeeze(mri.vol(xyz(elecId,1),:,:)),[mn mx]);
        axis square;
        hold on;
        
        if universalYes(anatOverlay)
            % Plot segmentation
            for a=1:sVol(2),
                for b=1:sVol(3),
                    if seg.vol(xyz(elecId,1),a,b)
                        segId=find(tbl{1}==seg.vol(xyz(elecId,1),a,b));
                        tempRgb=double([tbl{3}(segId) tbl{4}(segId) tbl{5}(segId)])/255;
                        hM=patch([-.5 .5 .5 -.5]+b,[-.5 -.5 .5 .5]+a,tempRgb);
                        set(hM,'LineStyle','none','FaceAlpha',0.3);
                    end
                end
            end
        end
        
        hm(2)=plot(xyz(elecId,3),xyz(elecId,2),'r.');
        set(hm(2),'color',elecRgb(elecId,:),'markersize',markerSize);
        %find image limits
        mxX=max(squeeze(mri.vol(xyz(elecId,1),:,:)),[],2);
        mxY=max(squeeze(mri.vol(xyz(elecId,1),:,:)),[],1);
        limXa=max(intersect(1:(sVol(3)/2),find(mxX==0)));
        limXb=min(intersect((sVol(3)/2:sVol(3)),find(mxX==0)));
        limYa=max(intersect(1:(sVol(2)/2),find(mxY==0)));
        limYb=min(intersect((sVol(2)/2:sVol(2)),find(mxY==0)));
        %keep image square
        tempMin=min([limXa limYa]);
        tempMax=max([limXb limYb]);
        if tempMin<tempMax,
            axis([tempMin tempMax tempMin tempMax]);
        end
        set(gca,'xtick',[],'ytick',[],'xdir','reverse');
        
        
        %subplot(133);
        axes('position',[xStart+wDelt*2 yStart wdth ht]);
        imagesc(squeeze(mri.vol(:,:,xyz(elecId,3))),[mn mx]);
        axis square;
        hold on;
        
        if universalYes(anatOverlay)
            % Plot segmentation
            for a=1:sVol(1),
                for b=1:sVol(2),
                    if seg.vol(a,b,xyz(elecId,3))
                        segId=find(tbl{1}==seg.vol(a,b,xyz(elecId,3)));
                        tempRgb=double([tbl{3}(segId) tbl{4}(segId) tbl{5}(segId)])/255;
                        hM=patch([-.5 .5 .5 -.5]+b,[-.5 -.5 .5 .5]+a,tempRgb);
                        set(hM,'LineStyle','none','FaceAlpha',0.3);
                    end
                end
            end
        end
        
        hm(3)=plot(xyz(elecId,2),xyz(elecId,1),'r.');
        set(hm(3),'color',elecRgb(elecId,:),'markersize',markerSize);
        %find image limits
        mxX=max(squeeze(mri.vol(:,:,xyz(elecId,3))),[],2);
        mxY=max(squeeze(mri.vol(:,:,xyz(elecId,3))),[],1);
        limXa=max(intersect(1:(sVol(3)/2),find(mxX==0)));
        limXb=min(intersect((sVol(3)/2:sVol(3)),find(mxX==0)));
        limYa=max(intersect(1:(sVol(2)/2),find(mxY==0)));
        limYb=min(intersect((sVol(2)/2:sVol(2)),find(mxY==0)));
        %keep image square
        tempMin=min([limXa limYa]);
        tempMax=max([limXb limYb]);
        if tempMin<tempMax,
            axis([tempMin tempMax tempMin tempMax]);
        end
        set(gca,'xtick',[],'ytick',[]);
        
        % Get anatomical label if aparc-file exists
        if exist(fullfile(mriDir,'aparc+aseg.mgz'),'file')
            anatLabel=vox2Seg(xyz(elecId,:),fsSub);
        else
            anatLabel = 'NA';
        end
        
        % Remove first 3 characters that indicate hemisphere and electrode
        % type
        formattedLabel=elecLabels{elecId}(4:end);
        formattedLabel=rmChar(formattedLabel,'_'); % remove underscore between electrode stem and #
        
        if universalYes(fullTitle)
            ht=textsc2014([formattedLabel '; mgrid coords(' num2str(elecMatrix(elecId,:)-1) '); fsurf coords(' num2str(xyz(elecId,:)) '); ' anatLabel], ...
                'title');
            set(ht,'fontsize',14,'fontweight','bold');
        else
            ht=textsc2014([formattedLabel '; Anatomical Location: ' anatLabel], ...
                'title');
            set(ht,'fontsize',16,'fontweight','bold');
        end
        set(ht,'position',[.5 .97 0]);
        
        if ~universalNo(printFigs)
            % Make sure PICS directory exists
            if ischar(printFigs)
                outPath=printFigs; % user specified directory
            else
                erPath=fullfile(fsdir,fsSub,'elec_recon');
                outPath=fullfile(erPath,'PICS');
            end
            if ~exist(outPath,'dir')
                dirSuccess=mkdir(outPath);
                if ~dirSuccess,
                    error('Could not create directory %s',dirSuccess);
                end
            end
            
            drawnow;
            figFname=fullfile(outPath,sprintf('%s_%sSlices',fsSub,elecLabels{elecId}));
            fprintf('Exporting figure to %s\n',figFname);
            print(figId,figFname,'-djpeg');
            pause(1);
            %close(figId);
        end
        
        if universalYes(pauseOn)
            fprintf('Paused. Press any key for next electrode.\n');
            pause;
        end
    end
end
fprintf('Done showing all electrodes.\n');