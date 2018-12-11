subj='PT001';
rootDir='/Users/davidgroppe/Desktop/';
task='tactile';
%globalFsDir
fsDir=getFsurfSubDir();
fsSubDir=fullfile(fsDir,subj);
elecReconDir=fullfile(fsSubDir,'elec_recon');
sessionId=1;
% Get Freesurfer Version
fsurfVersionFname=fullfile(fsSubDir,'scripts','build-stamp.txt');
fid=fopen(fsurfVersionFname);
fsurfVersion=fgetl(fid);
fclose(fid);

% BIDS directories
taskDir=fullfile(rootDir,sprintf('%s_task',task));
[SUCCESS,MESSAGE,MESSAGEID] = mkdir(taskDir);
bidsSubDir=fullfile(taskDir,sprintf('sub-%s',subj));
[SUCCESS,MESSAGE,MESSAGEID] = mkdir(bidsSubDir);
ieegDir=fullfile(bidsSubDir,'ieeg');
[SUCCESS,MESSAGE,MESSAGEID] = mkdir(ieegDir);
derivDir=fullfile(taskDir,'derivatives');
[SUCCESS,MESSAGE,MESSAGEID] = mkdir(derivDir);
derivSubDir=fullfile(derivDir,sprintf('sub-%s',subj));
[SUCCESS,MESSAGE,MESSAGEID] = mkdir(derivSubDir);


%% TODO change CT to POSTIMP


%% TODO map to avg brain


%% Export neuroimaging to anat folder
anatDir=fullfile(bidsSubDir,'anat');
[SUCCESS,MESSAGE,MESSAGEID] = mkdir(anatDir);
copyfile(fullfile(elecReconDir,'T1.nii.gz'),fullfile(anatDir,'T1_fsurf.nii.gz'));
postimpRawFname='postimpRaw.nii.gz';
if ~exist(postimpRawFname,'file')
    postimpRawFname='postopCT.nii.gz'; % This was the original filename for the postimplant scan (CT or MRI)
end
copyfile(fullfile(elecReconDir,postimpRawFname),fullfile(anatDir,'postimp_raw.nii.gz'));


%% Export neuroimaging to derivatives folder
% copy FreeSurfer version to derivatives folder
copyfile(fsurfVersionFname,fullfile(derivSubDir,'freesurferVersion.txt'));

% copy postimplant scan aligned to preimplant scan
postimpAlignedFname='postInPre.nii.gz';
if ~exist(postimpAlignedFname,'file')
    postimpAlignedFname='ctINt1.nii.gz'; % This was the original filename for the postimplant scan (CT or MRI)
end
copyfile(fullfile(elecReconDir,postimpAlignedFname),fullfile(derivSubDir,'postInPre.nii.gz'));

% Copy anatomical labels
labelFiles={'aparc.a2009s.annot','aparc.annot'};
tempSubDir='label';
fsurf2bids(labelFiles,fullfile(fsSubDir,tempSubDir),fullfile(derivSubDir,tempSubDir));

% Copy surface files to BIDS derivatives dir
% DG: I don't know if all of these files are necessary
surfFiles={'area.pial','curv','curv.pial','inflated','pial','pial-outer-smoothed','sphere','sphere.reg','white'};
tempSubDir='surf';
fsurf2bids(surfFiles,fullfile(fsSubDir,tempSubDir),fullfile(derivSubDir,tempSubDir));

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

% TODO Yeo atlas if the files exist


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


%% Import electrode locations
%TODO make this able to deal with "ct" or "postimplant" fnames
coordTypes={'ct','inf','leptp','leptovox','pial','pialvox'};
nCoordSpaces=length(coordTypes);
for sLoop=1:1, %% ?? debuging
    %for sLoop=1:nCoordSpaces,
    ielvisPostfix=upper(coordTypes{sLoop});
    
    % Import coordinates
    xyzFname=fullfile(elecReconDir,sprintf('%s.%s',subj,ielvisPostfix));
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
    splitHdr=strsplit(firstLine,9); % split on tabs
    if length(splitHdr)>=2,
        % method is specified in header
        brainShiftCorrectMethod=splitHdr{2}; 
    else
        % figure out which brain shift correction method was used based on log
        splitHdr2=strsplit(firstLine,32);
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
    elecTsvFname=fullfile(ieegDir,sprintf('sub-%s_ses-%.2d_space-%s_electrodes.tsv', ...
        subj,sessionId,elecSpace));
    fid=fopen(elecTsvFname,'w');
    % Header
    fprintf(fid,'name\tx\ty\tz\tsize\ttype\tmanufacturer\n');
    % Elec info
    for eLoop=1:nElec,
        elecManufacturer='n/a'; % TODO make it possible to import this from file
        elecSize=nan; % TODO make it possible to import this from file
        % TODO add hemisphere to this?
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
        fprintf(fid,'%d\t%f\t%f\t%f\t%f\t%s\t%s\n',elecNames{eLoop},xyzCoord(eLoop,1), ...
            xyzCoord(eLoop,2),xyzCoord(eLoop,3),elecSize,thisElecType,elecManufacturer);
    end
    fclose(fid);
    
    
    %% Create coordsystem.json file
    coordJsonFname=fullfile(ieegDir,sprintf('sub-%s_ses-%.2d_space-%s_coordsystem.json', ...
        subj,sessionId,elecSpace));
    anatFile='??'; % TODO fix
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
    fprintf(fid,'"IntendedFor": "%s"\n',anatFile);
    fprintf(fid,'}\n');
    fclose(fid);
    
end


