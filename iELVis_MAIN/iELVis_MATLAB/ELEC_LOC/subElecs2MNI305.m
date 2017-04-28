function [avgCoords, elecNames]=subElecs2MNI305(subj,createTxtFile)
%function [avgCoords, elecNames]=subElecs2MNI305(subj,createTxtFile)
%
% This function takes electrode RAS LEPTO coordinates (i.e., patient-specific
% brain shift corrected coordinates) and maps them to the corresponding 
% location in MNI305 space using an AFFINE transformation. As described as
% follows on the FreeSurfer wiki:
%
% http://surfer.nmr.mgh.harvard.edu/fswiki/CoordinateSystems
% "2. I have an RAS point on the surface (tkrR tkrA tkrS) ("Vertex RAS" from
% tksurfer) and want to compute the MNI305 RAS that corresponds to this point:
% 
% MNI305RAS = TalXFM*Norig*inv(Torig)*[tkrR tkrA tkrS 1]'
% TalXFM: subject/orig/transforms/talairach.xfm Norig: mri_info --vox2ras orig.mgz Torig: mri_info --vox2ras-tkr orig.mgz 
% 
% Test: click on a point in tksurfer. Use "Vertex RAS" to compute MNI305RAS. 
% Compare MNI305RAS to "Vertex MNI Talairach". Also, hit Save Point in tksurfer. 
% In tkmedit, Goto Saved Point, compare MNI305RAS to "MNI Coordinates".
% Note: you may need to use the "Volume RAS" from tkmedit in the computation 
% above to get an exact match."
%
%
% Required Input:
%   subj = FreeSurfer subject name
%
% Optional Input:
%   createTxtFile = [1 or 0] If 1, the MNI305 coordinates are saved to a
%                   *.MNI305 text file in the "elec_recon" subfolder of the 
%                   patient's FreeSurfer subject directory.
%
% Outputs:
%   avgCoords = Electrode coordinates in MNI305 space that can be visualized
%               on the FreeSurfer avg brain (RAS coordinates)
%   elecNames = Channel names
%
%
% Authors:
% Pierre Megevand & David Groppe
% Mehtalab, Feinstein Institute for Medical Research
% Honeylab, University of Toronto
% 2015-2016
%

% TODO
% -check output coordinates
% -make creating MNI305 file optional


%% Input Arguments
if nargin<2,
    createTxtFile=0;
end

%% FreeSurfer Subject Directory
fsDir=getFsurfSubDir();
subDir=fullfile(fsDir,subj);

if ~exist(subDir,'dir')
    error('Folder for %s is not present in FreeSurfer subjects directory (%s).',subj,fsDir);
end


% load MRI header
mgzFname=fullfile(subDir,'mri','orig.mgz');
checkFile(mgzFname);
MRIhdr=MRIread(mgzFname,true);
Norig=MRIhdr.vox2ras; 
Torig=MRIhdr.tkrvox2ras;


%%
TalTransform=freesurfer_read_talxfm(fullfile(subDir,'mri','transforms','talairach.xfm'));


%% Import Electrode Coordinates in Patient Space and Electrode Names
ptntCoordFile=fullfile(subDir,'elec_recon',[subj '.LEPTO']);
tempCsv=csv2Cell(ptntCoordFile,' ',2);
nElec=size(tempCsv,1);
ptntCoords=zeros(nElec,3);
for a=1:nElec,
   for b=1:3,
    ptntCoords(a,b)=str2double(tempCsv{a,b});
   end
end

% Get electrode names
elecNamesFile=fullfile(subDir,'elec_recon',[subj '.electrodeNames']);
tempCsv=csv2Cell(elecNamesFile,' ',2);
elecNames=tempCsv(:,1);


%%
% Pierre took this code from http://surfer.nmr.mgh.harvard.edu/fswiki/CoordinateSystems
% Section #2 of "Transforms within a subject's anatomical space"
% "2. I have an RAS point on the surface (tkrR tkrA tkrS) ("Vertex RAS"
% from tksurfer) and want to compute the MNI305 RAS that corresponds to this point:"
avgCoords=(TalTransform*Norig*(Torig\[ptntCoords'; ones(1, nElec)]))';
avgCoords=avgCoords(:,1:3);

%% Get LEPTOVOX file time stamp
leptoFileId=fopen(ptntCoordFile,'r');
timeStamp=fgetl(leptoFileId);
fclose(leptoFileId);

%% Write coordinates to file
if universalYes(createTxtFile),
    outFname=fullfile(subDir,'elec_recon',[subj '.MNI305']);
    fprintf('Saving MNI305 coordinates to file %s\n',outFname);
    fid=fopen(outFname,'w');
    fprintf(fid,'%s\n',timeStamp);
    fprintf(fid,'R A S\n');
    for a=1:nElec,
        fprintf(fid,'%f %f %f\n',avgCoords(a,1),avgCoords(a,2),avgCoords(a,3));
    end
    fclose(fid);
end
