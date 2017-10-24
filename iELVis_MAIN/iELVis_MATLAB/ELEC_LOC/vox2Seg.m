function [anatLabel, ROIs]=vox2Seg(coordILA,fsSub,atlas,prob)
%function anatLabel=vox2Seg(coordILA,fsSub)
%
% Inputs:
%  coordILA - 3D vector indicating the coordinates of a voxel in a
%             FreeSurfer MRI. First coordinate is Sup->Inf (i.e., 1=most
%             superior slice). Second coordinate is R->L (i.e., 1=the
%             rightmost slice). Third coordinate is P->A (i.e., 1=the
%             most posterior slice).
%  fsSub    - The FreeSurfer folder for the patient
%
%  atlas    - Anatomical atlas to use:
%           'DK'=Desikan-Killiany {default} (from aparc+aseg.mgz)
%           'D' =Destrieux                  (from aparc.a2009+aseg.mgz)
%
%  prob     - Binary (1/'yes' or 0/'no'
%            If non-zero, anatLabel will be run for voxels in the vicinity
%            of the input coordinates to estimate the representativeness of
%            the anatomical label (offset = 3x3 voxels)
%
% Output:
%  anatLabel - The part of the brain that voxel belongs to (e.g., Left
%              Hippocampus) according to aparc+aseg.mgz
%
%  COIs      - three fields: anatLabel and corresponding prob in vicinity
%                            name for the anatomical label in the vicinity
%                            count give the proportion of voxel labelled with each name
%
%  Example:
%  >>coordILA=[147 144 115];
%  >>anatLabel=vox2Seg(coordILA,'NiAs')
%
% Author: David M. Groppe
% Feb. 2015
% Feinstein Institute for Medical Research/Univ. of Toronto
%
% Sept. 2017 - Manuel R. Mercier (manuel.mercier@a3.epfl.ch) from CerCo lab (CNRS)
% modifs:
% Add atlas option to choose between atlas  uses the Destrieux atlas instead of Desikan-Killiany. 
%

% default atlas
if nargin<3
    atlas='DK';
end

offset = 2;
if nargin<4
    prob='n';
end

% Load aseg volume
fsdir=getFsurfSubDir();

switch upper(atlas)
    case 'DK'
        asegFname=fullfile(fsdir, fsSub, '/mri/aparc+aseg.mgz');
    case 'D'
        asegFname=fullfile(fsdir, fsSub, '/mri/aparc.a2009s+aseg.mgz');
    otherwise
        warning('only Desikan-Killiany and Destrieux atlases are available');
end
%asegFname=[fsdir '/' fsSub '/mri/aseg.mgz'];
<<<<<<< HEAD

=======
asegFname=fullfile(fsdir,fsSub,'mri','aparc+aseg.mgz');
%asegFname=[fsdir '/' fsSub '/mri/aparc.a2009s+aseg.mgz'];
>>>>>>> refs/remotes/iELVis/master
if ~exist(asegFname,'file')
   error('File %s not found.',asegFname); 
end

aseg=MRIread(asegFname);

%% Load table
pathstr = fileparts(which('mgrid2matlab'));
inFile=fullfile(pathstr,'FreeSurferColorLUTnoFormat.txt');
if ~exist(inFile,'file')
    error('Could not find file %s',inFile);
end
fid=fopen(inFile,'r');
%fid=fopen('/Applications/freesurfer/FreeSurferColorLUTnoFormat.txt','r');
tbl=textscan(fid,'%d%s%d%d%d%d');
fclose(fid);

%% Find anatomical region corresponding to voxel
id=find(aseg.vol(coordILA(1),coordILA(2),coordILA(3))==tbl{1});
% id=min(id);
anatLabel=tbl{2}{id};

%% Get anatomical labels for voxel in the vicinity and compute their recurrence 
if universalYes(prob)
    % get the euclidean distances between the electrode and every voxel in the MRI
    for i=1:aseg.volsize(1)
        for j=1:aseg.volsize(2)
            for k=1:aseg.volsize(3)
                distances(i,j,k)=sqrt((i-coordILA(1))^2+(j-coordILA(2))^2+(k-coordILA(3))^2); % could be replaced by ~ norm([i j k] - [x y z])
            end
        end
    end
    % do not include the offset
    tmp=aseg.vol(distances<offset);
    % find the regions
    tmp_region ={};
    for i = 1:length(tmp)
        tmp_region{i,1} =tbl{2}{tbl{1}==tmp(i)};
    end
    ROIs.name  = unique(tmp_region);
    for i = 1:length(ROIs.name)
        ROIs.count(i) = sum(count(tmp_region,ROIs.name{i}));
    end
    mainLab = find(contains(ROIs.name,anatLabel));
    ROIs.center{2} = num2str(ROIs.count(mainLab)/sum(ROIs.count));
    ROIs.center{1} = anatLabel;
    ROIs.offset    = offset;
else
    ROIs=[];
end
