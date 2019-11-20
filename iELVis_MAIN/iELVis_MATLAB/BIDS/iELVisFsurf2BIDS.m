function iELVisFsurf2BIDS(subj,bidsRootDir,sessionId, atlas, elecReconDir, exportAnat, exportDeriv)
% function iELVisFsurf2BIDS(subj,bidsRootDir,sessionId, atlas, elecReconDir, exportAnat, exportDeriv)
% 
% Takes iELVis electrode location and neuroimaging information and copies
% it into a new directory according to BIDS-iEEG conventions.
%
% Here is a preprint of the paper introducing BIDS-iEEG and describing its 
% format:
%     https://psyarxiv.com/r7vc2/
%
% Required Inputs:
%  subj - [string] Patient FreeSurfer ID
%
%  bidsRootDir - [string] The root BIDS directory into which all data from a study
%                will be stored
%
% Optional Inputs:
%  sessionId - [integer or string] A string or number indicating the "session" 
%              in which the corresponding iEEG data have been recorded. This 
%              will be used to create filenames of electrode locations as it is
%              possible that data may be recorded in multiple sessions and
%              the electrodes recorded from may differ across
%              sessions. If left as empty (i.e. [] ) then no session folder
%              will be created and no session name will be added to
%              electrode file names
%              {default: []}
%
%  atlas - Anatomical atlas to use for identifying the anatomical region
%          nearest each electrode. Options:
%             'DK'=Desikan-Killiany {default}
%             'D' =Destrieux
%             'Y7'=Yeo 7-network resting state fMRI atlas
%             'Y17'=Yeo 17-network resting state fMRI atlas
%           Left empty (i.e. [] ) will use default 'DK' atlas
%
%  elecReconDir - [string] The directory from which iELVis information
%                (e.g., electrode names and coordinates) will be taken 
%                from. Leaving empty (i.e. [] ) will use default 
%                {default: the "elec_recon" subfolder of the patient's
%                FreeSurfer subject dir}
%               
%  exportAnat - [boolean] Whether a sub/anat directory is created and
%               populated with pre and post implant MRI scans
%               {default: 1}
%
%  exportDeriv - [boolean] Whether a derivatives/sub subdirectory is
%               created and populated with freesurfer-derived files that
%               contributed to reconstruction
%               {default: 1}
%
% Examples:
%  iELVisFsurf2BIDS('PT001','~/Desktop/HandMotor');
%       - Electrode file names will be 'sub-PT001_space-....'
%
%  iELVisFsurf2BIDS('PT001','~/Desktop/HandMotor', 'implant', 'D', [], 0, 0);
%       - Electrode file names will be 'sub-PT001_ses-implant_space-...'
%       - Anatomic locations are based on Destrieux atlas
%       - Uses default elec_recon directory
%       - Does not create 'derivatives/iElVis/sub-PT001' or
%           'sub-PT001/ses-implant/anat' directories
% 

%extraElecInfoFname=[]; % TODO make this an argument to the function

fprintf('Exporting data from %s into iEEG-BIDS format.\n',subj);
fprintf('iEEG-BIDS root directory is %s\n',bidsRootDir);

% Get FreeSurfer directories
fsDir=getFsurfSubDir();
fsSubDir=fullfile(fsDir,subj);

if nargin<3 || isempty(sessionId),
    sessionId=[];
    filesBaseName = sprintf('sub-%s', subj);
else
   if isnumeric(sessionId)
      sessionId = sprintf('%02d', sessionId); 
   end
   filesBaseName = sprintf('sub-%s_ses-%s', subj,sessionId);
end
if nargin<4,
    atlas=[];
end
if nargin<5 || isempty(elecReconDir),
    elecReconDir=fullfile(fsSubDir,'elec_recon');
end
if nargin < 6,
    exportAnat = 1;
end
if nargin < 7,
    exportDeriv = 1;
end
%% Get Freesurfer Version
fsurfVersionFname=fullfile(fsSubDir,'scripts','build-stamp.txt');
if exist(fsurfVersionFname,'file')
    fid=fopen(fsurfVersionFname);
    fsurfVersion=fgetl(fid);
    fclose(fid);
else
    fsurfVersion='FreeSurferVersionUnknown';
end


%% Define/Create BIDS directories
if ~isfolder(bidsRootDir); [SUCCESS,MESSAGE,MESSAGEID] = mkdir(bidsRootDir); end

if ~isempty(sessionId)
    bidsSubDir=fullfile(bidsRootDir,sprintf('sub-%s%sses-%s',subj, filesep, sessionId));
else
    bidsSubDir=fullfile(bidsRootDir,sprintf('sub-%s',subj));
end

if ~isfolder(bidsSubDir); [SUCCESS,MESSAGE,MESSAGEID] = mkdir(bidsSubDir); end
ieegDir=fullfile(bidsSubDir,'ieeg');
if ~isfolder(ieegDir); [SUCCESS,MESSAGE,MESSAGEID] = mkdir(ieegDir); end
derivDir=fullfile(bidsRootDir,'derivatives','iELVis');
if ~isfolder(derivDir) && istrue(exportDeriv); [SUCCESS,MESSAGE,MESSAGEID] = mkdir(derivDir); end
derivSubDir=fullfile(derivDir,sprintf('sub-%s',subj));
if ~isfolder(derivSubDir) && istrue(exportDeriv); [SUCCESS,MESSAGE,MESSAGEID] = mkdir(derivSubDir); end


%% Export neuroimaging to anat folder
if istrue(exportAnat)
    anatDir=fullfile(bidsSubDir,'anat');
    [SUCCESS,MESSAGE,MESSAGEID] = mkdir(anatDir);
    copyfile(fullfile(fsSubDir,'mri','orig','001.mgz'),fullfile(anatDir,'preimpRaw.mgz'));
    postimpRawFname='postimpRaw.nii.gz';
    if ~exist(fullfile(elecReconDir,postimpRawFname),'file')
        % Data must have been processed with original version of iELVis that
        % did not standardize the raw postimplant nii filename.
        % Ask user to select the postimplant nii file
        homeDir=pwd;
        cd(elecReconDir);
        disp('Select postimplant raw nii file');
        [postimpRawFname, tempPathName] = uigetfile('*.nii.gz;*.nii', 'Select postimplant raw nii file');
        if postimpRawFname==0
           error('You need to select a file that contains the raw postimplant volume.'); 
        end
        cd(homeDir);
        copyfile(fullfile(tempPathName,postimpRawFname),fullfile(anatDir,'postimpRaw.nii.gz'));
    else
        copyfile(fullfile(elecReconDir,postimpRawFname),fullfile(anatDir,'postimpRaw.nii.gz'));
    end
end

%% Import information about electrode size and manufacturer TODO make this work once iEEG-BIDS format is set
% if ~isempty(extraElecInfoFname)
%     extraElecInfo=csv2Cell(extraElecInfoFname,',',1);
%     extraElecNames=extraElecInfo(:,1);
%     extraElecManufacturer=extraElecInfo(:,1);
%     nExtra=size(extraElecNames,1); 
% end


%% Export neuroimaging to derivatives folder
% copy FreeSurfer version to derivatives folder
if istrue(exportDeriv)
    copyfile(fsurfVersionFname,fullfile(derivSubDir,'freesurferVersion.txt'));

    % Copy FreeSurfer pre-processed preimplant scan
    preimpPreprocFname='T1.nii.gz';
    copyfile(fullfile(elecReconDir,preimpPreprocFname),fullfile(derivSubDir,'preimplantFsurf.nii.gz'));

    % Copy postimplant scan aligned to preimplant scan
    postimpAlignedFname='postInPre.nii.gz';
    if ~exist(postimpAlignedFname,'file')
        postimpAlignedFname='ctINt1.nii.gz'; % This was the original filename for the postimplant scan (CT or MRI)
    end
    copyfile(fullfile(elecReconDir,postimpAlignedFname),fullfile(derivSubDir,'postInPre.nii.gz'));

    % Copy anatomical labels
    labelFiles={'aparc.a2009s.annot','aparc.annot'};
    tempSubDir='label';
    fsurf2BIDS(labelFiles,fullfile(fsSubDir,tempSubDir),fullfile(derivSubDir,tempSubDir));

    % Copy Yeo anatomical labels if they exist
    labelFiles={'Yeo2011_17Networks_N1000.mat','Yeo2011_7Networks_N1000.mat'};
    tempSubDir='label';
    fsurf2BIDS(labelFiles,fullfile(fsSubDir,tempSubDir),fullfile(derivSubDir,tempSubDir),1);

    % Copy surface files to BIDS derivatives dir
    % DG: I don't know if all of these files are necessary
    surfFiles={'area.pial','curv','curv.pial','inflated','pial','pial-outer-smoothed','sphere','sphere.reg','white'};
    tempSubDir='surf';
    fsurf2BIDS(surfFiles,fullfile(fsSubDir,tempSubDir),fullfile(derivSubDir,tempSubDir));

    % Copy mri volumes to BIDS derivatives dir
    mriFiles={'aparc+aseg.mgz','brainmask.mgz','orig.mgz'};
    tempSubDir='mri';
    derivMriDir=fullfile(derivSubDir,tempSubDir);
    [SUCCESS,MESSAGE,MESSAGEID] = mkdir(derivMriDir);
    for a=1:length(mriFiles),
        copyfile(fullfile(fsSubDir,tempSubDir,mriFiles{a}),fullfile(derivMriDir,mriFiles{a}));
    end
    tempSubDir='transforms';
    derivTransDir=fullfile(derivMriDir,tempSubDir);
    [SUCCESS,MESSAGE,MESSAGEID] = mkdir(derivTransDir);
    fname='talairach.xfm';
    copyfile(fullfile(fsSubDir,'mri',tempSubDir,fname),fullfile(derivTransDir,fname));

end
%% Export FreeSurfer avg brain and associate info to derivates
if istrue(exportDeriv)
    
    derivFsavgDir=fullfile(derivDir,'sub-fsaverage');
    [SUCCESS,MESSAGE,MESSAGEID] = mkdir(derivFsavgDir);
    derivFsavgSurfDir=fullfile(derivFsavgDir,'surf');
    [SUCCESS,MESSAGE,MESSAGEID] = mkdir(derivFsavgSurfDir);

    fsurfAvgSurfDir=fullfile(fsDir,'fsaverage','surf');
    fnameStems={'inflated','pial'};
    for a=1:2,
        if a==1
            hem='lh';
        else
            hem='rh';
        end
        for b=1:length(fnameStems),
            copyfile(fullfile(fsurfAvgSurfDir,[hem '.' fnameStems{b}]), ...
                fullfile(derivFsavgSurfDir,[hem '.' fnameStems{b}]));
        end
    end

    derivFsavgLabelDir=fullfile(derivFsavgDir,'label');
    [SUCCESS,MESSAGE,MESSAGEID] = mkdir(derivFsavgLabelDir);
    fsurfAvgLabelDir=fullfile(fsDir,'fsaverage','label');
    fnameStems={'aparc.annot','aparc.a2009s.annot','Yeo2011_17Networks_N1000.annot', ...
        'Yeo2011_7Networks_N1000.annot'};
    for a=1:2,
        if a==1
            hem='lh';
        else
            hem='rh';
        end
        for b=1:length(fnameStems),
            copyfile(fullfile(fsurfAvgLabelDir,[hem '.' fnameStems{b}]), ...
                fullfile(derivFsavgLabelDir,[hem '.' fnameStems{b}]));
        end
    end
end  

%% Import electrode names, type, and hemisphere
% Electrode names
elecNamesFname=fullfile(elecReconDir,sprintf('%s.electrodeNames',subj));
elecNamesCsv=csv2Cell(elecNamesFname,' ',2);
nElec=size(elecNamesCsv,1);
elecNames=cell(nElec,1);
elecType=cell(nElec,1);
elecHem=cell(nElec,1);
for a=1:nElec,
    elecNames{a}=elecNamesCsv{a,1};
    elecType{a}=elecNamesCsv{a,2};
    elecHem{a}=elecNamesCsv{a,3};
end


%% Map to avg brain if it hasn't been done already
avgCoordFname=fullfile(elecReconDir,[subj '.FSAVERAGE']);
if ~exist(avgCoordFname,'file')
    fprintf('Creating text file of electrodes coordinates in FreeSurfer avg. brain space.\n');
    cfg=[];
    cfg.plotEm=0;
    sub2AvgBrain(subj,cfg);
end;


%% Get anatomical labels for electrodes if they don't already exist
elecAnatLabel=cell(nElec,1);
if isempty(atlas)
    for a=1:nElec,
        elecAnatLabel{a}='n/a';
    end
else
    atlas=upper(atlas);
    anatLabelsFname=fullfile(elecReconDir,[subj '_' atlas '_AtlasLabels.tsv']);
    if ~exist(anatLabelsFname,'file')
        fprintf('Creating text file that indicates which anatomical region each electrode is closest to.\n');
        fprintf('%s\n',anatLabelsFname);
        elec2Parc(subj,atlas,1);
    end
    anatLabelsTable=csv2Cell(anatLabelsFname,9);
    for a=1:nElec,
        tempId=findStrInCell(elecNames{a},anatLabelsTable(:,1),1);
        elecAnatLabel{a}=anatLabelsTable{tempId,2};
    end
end

%% Import electrode locations and create channels.tsv and coordsystem.json file
% Note this can deal with "ct" or "postimplant" fnames
coordTypes={'postimplant','inf','lepto','leptovox','pial','pialvox','fsaverage'};
nCoordSpaces=length(coordTypes);
for sLoop=1:nCoordSpaces,
    ielvisPostfix=upper(coordTypes{sLoop});
    
    % Import coordinates
    xyzFname=fullfile(elecReconDir,sprintf('%s.%s',subj,ielvisPostfix));
    if strcmpi(ielvisPostfix,'POSTIMPLANT')
        if ~exist(xyzFname,'file')
            % Try older filename convention
            xyzFname=fullfile(elecReconDir,sprintf('%s.CT',subj));
            if ~exist(xyzFname,'file')
                error('Missing *.POSTIMPLANT/*.CT file for subject %s\n',subj);
            end
        end
    end
    xyzCoordStr=csv2Cell(xyzFname,' ',2);
    if size(xyzCoordStr,1)~=nElec,
        error('# of electrodes in %s does not match that in %s\n',xyzFname,elecNamesFname);
    end
    % Convert coordinates from string to double
    xyzCoord=zeros(nElec,3);
    for a=1:nElec,
        for b=1:3,
            xyzCoord(a,b)=str2double(xyzCoordStr{a,b});
        end
    end
    
    % Import brain shift correction method
    fid=fopen(xyzFname,'r');
    firstLine=fgetl(fid);
    fclose(fid);
    splitHdr=strsplit(firstLine,'\t'); % split on tabs
    if length(splitHdr)>=2,
        % method is specified in header
        brainShiftCorrectMethod=splitHdr{2}; 
    else
        % figure out which brain shift correction method was used based on log
        splitHdr2=strsplit(firstLine,'\t');
        dateGenerated=datetime(splitHdr2{1});
        logDate=datestr(dateGenerated,'yyyy-mm-dd');
        logFname=sprintf('localization_process_%s.log',logDate);
        logFname=fullfile(elecReconDir,logFname);
        if exist(logFname,'file')
            fid=fopen(logFname,'r');
            tempLine=fgetl(fid);
            tempLine=fgetl(fid);
            tempLine=fgetl(fid);
            if strfind('Dykstra',tempLine)
                brainShiftCorrectMethod='dykstra-preIeegBids';
            else
                brainShiftCorrectMethod='yangWang-preIeegBids';
            end
            fclose(fid);
        else
            brainShiftCorrectMethod='BrainShiftCorrectionMethodUnknown';
        end
    end
    
    %% Create channels.tsv file
    elecSpace=coordTypes{sLoop};
    % outFname=fullfile(ieegDir,sprintf('sub-%.2d_task-%s_run-%.2_channels.tsv', ...
    %     subId,task,runId));
    elecTsvFname=fullfile(ieegDir,sprintf('%s_space-%s_electrodes.tsv', ...
        filesBaseName,elecSpace));
    fid=fopen(elecTsvFname,'w');
    % Header
    fprintf(fid,'name\tx\ty\tz\themisphere\tsize\ttype\tmanufacturer\tanatomicalLocation\n');
    % Elec info
    for eLoop=1:nElec,
        elecManufacturer='n/a'; % TODO make it possible to import this from file
        elecSize=nan; % TODO make it possible to import this from file
        switch elecType{eLoop}
            case 'D'
                thisElecType='depth';
            case 'G'
                thisElecType='grid';
            case 'S'
                thisElecType='strip';
            otherwise
                error('Unrecognized value of elecType: %s',elecType{eLoop});
        end
        fprintf(fid,'%s\t%f\t%f\t%f\t%s\t%f\t%s\t%s\t%s\n',elecNames{eLoop},xyzCoord(eLoop,1), ...
            xyzCoord(eLoop,2),xyzCoord(eLoop,3),elecHem{eLoop},elecSize,thisElecType,elecManufacturer,elecAnatLabel{eLoop});
        %         fprintf(fid,'%s\t%f\t%f\t%f\t%f\t%s\t%s\n',elecNames{eLoop},xyzCoord(eLoop,1), ...
        %             xyzCoord(eLoop,2),xyzCoord(eLoop,3),elecSize,thisElecType,elecManufacturer);
    end
    fclose(fid);
    
    
    %% Create coordsystem.json file
    coordJsonFname=fullfile(ieegDir,sprintf('%s_space-%s_coordsystem.json', ...
        filesBaseName,elecSpace));
    switch ielvisPostfix
        case 'FSAVERAGE'
            anatFile=sprintf('/derivatives/iELVis/fsaverage/surf/*.pial');
        case {'LEPTOVOX','PIALVOX'}
            anatFile=sprintf('/sub-%s/anat/T1_fsurf.nii.gz',subj);
        case 'INF'
            anatFile=sprintf('/derivatives/sub-%s/surf/*.inflated',subj);
        otherwise
            anatFile=sprintf('/derivatives/sub-%s/surf/*.pial',subj);
    end
    fid=fopen(coordJsonFname,'w');
    fprintf(fid,'{\n');
    fprintf(fid,'"iEEGCoordinateSystem": "%s",\n',elecSpace);
    fprintf(fid,'"iEEGCoordinateUnits": "mm",\n');
    fprintf(fid,'"iEEGCoordinateProcessingDescription": "%s",\n',brainShiftCorrectMethod);
    if strfind('dykstra',brainShiftCorrectMethod),
        brainShiftReference='Dykstra et al., 2011 NeuroImage; Groppe et al., 2017 JNeuroMeth';
    else
       brainShiftReference='Yang, Wang et al., 2012 NeuroImage; Groppe et al., 2017 JNeuroMeth';
    end
    fprintf(fid,'"iEEGCoordinateProcessingReference": "%s",\n',brainShiftReference);
    fprintf(fid,'"t1ProcessingDescription": "%s",\n',fsurfVersion);
    fprintf(fid,'"t1ProcessingReference": "http://freesurfer.net/fswiki/FreeSurferMethodsCitation",\n');
    fprintf(fid,'"IntendedFor": "%s",\n',anatFile);
    if ~isempty(atlas)
        switch atlas
            case 'DK'
                atlasLong='Desikan-Killiany';
            case 'D'
                atlasLong='Destrieux';
            case 'Y7'
                atlasLong='Yeo 7-network';
            case 'Y17'
                atlasLong='Yeo 17-network';
        end
        fprintf(fid,'"Atlas used for electrode anatomical labels": "%s"\n',atlasLong);
    end
    fprintf(fid,'}\n');
    fclose(fid);
end


