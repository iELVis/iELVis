function [avgCoords, elecNames, isLeft]=depths2AvgBrain(subj)
%function [avgCoords, elecNames, isLeft]=depths2AvgBrain(subj)
%
% This function takes RAS CT coordinates of depth electrodes (i.e.,
% the depth electrode coordinates in the postimplant CT or MRI) and maps it 
% to the corresponding location in MNI305 space using an affine
% transformation. These coordinates can then be used to visualize electrode
% locations on the FreeSurfer average brain.
%
% Inputs:
%   subj = FreeSurfer subject name
%
% Outputs:
%   avgCoords = Electrode coordinates in MNI305 space that can be visualized
%               on the FreeSurfer avg brain (RAS coordinates)
%   elecNames = Channel names with the participant's name appended to the
%               beginning (e.g., PT001-Gd1)
%   isLeft    = N-dimensional binary vector where N is the # of electrodes.
%               0 indicates that an electrode is on/in the right hemisphere.
%               1 indicates a left hemisphere location.
%
%
% Authors:
% Pierre Megevand & David Groppe
% Mehtalab, Feinstein Institute for Medical Research
% Honeylab, University of Toronto
% 2015-2016
%

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
ptntCoordFile=fullfile(subDir,'elec_recon',[subj '.CT']);
tempCsv=csv2Cell(ptntCoordFile,' ',2);
nElec=size(tempCsv,1);
ptntCoords=zeros(nElec,3);
for a=1:nElec,
   for b=1:3,
    ptntCoords(a,b)=str2double(tempCsv{a,b});
   end
end

elecNamesFile=fullfile(subDir,'elec_recon',[subj '.electrodeNames']);
tempCsv=csv2Cell(elecNamesFile,' ',2);
elecNames=tempCsv(:,1);
isLeft=zeros(nElec,1);
isDepth=zeros(nElec,1);
for a=1:nElec,
    if tempCsv{a,3}=='L'
        isLeft(a)=1;
    end
    if tempCsv{a,2}=='D'
        isDepth(a)=1;
    end
    %elecNames{a}=[subj '-' elecNames{a}];
end
% Select only depths
elecNames=elecNames(isDepth==1);
isLeft=isLeft(isDepth==1);
ptntCoords=ptntCoords(isDepth==1,:);
nElec=sum(isDepth);

%%
% Pierre took this code from http://surfer.nmr.mgh.harvard.edu/fswiki/CoordinateSystems
% Section #2 of "Transforms within a subject's anatomical space"
% "2. I have an RAS point on the surface (tkrR tkrA tkrS) ("Vertex RAS"
% from tksurfer) and want to compute the MNI305 RAS that corresponds to this point:"
avgCoords=(TalTransform*Norig*(Torig\[ptntCoords'; ones(1, nElec)]))';
avgCoords=avgCoords(:,1:3);
