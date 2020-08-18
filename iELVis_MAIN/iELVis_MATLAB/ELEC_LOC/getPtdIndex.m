function PTD_idx = getPtdIndex(fs_subj)
% function PTD_idx = getPtdIndex(fs_subj)
%
% Compute Proximal Tissue Density (PTD) for each electrode as described in
% Mercier et al., Neuroimage 2017
%
% PTD is an index reflecting the density of neocortical gray and white matter surrounding
% a stereotactic electrode that has its centroid either in the neocortical gray or in the
% white matter. It is: 
%
% PTD=(nGray-nWhite)/(nGray+nWhite)
% Where nGray is the # of neocortical gray matter voxels and nWhite is the
% number of white matter voxels in the 3x3 voxel cube centered on the contact. 
%
% If the contact centroid is not in white or neocortical gray matter, PTD =
% NaN.
%
% Be careful when a contact is in the vicinity of subcortical structures
% PTD is computed exclusively by taking into account surrounding neocortical
% gray and white matter voxels; no other tissue will be taken into account (e.g.
% Hippocampus, Amygdala...). Voxel labels are taken from the FreeSurfer
% aparc+aseg.mgz file.
%
% Please note that this function belongs to the iELVis toolbox
% and is therefore subjected to the same regulations
% (contacts: david.m.groppe@gmail.com ; manuel.mercier@a3.epfl.ch)
%
% input: fs_subj - name of the patient's FreeSurfer subject folder
%
% output: PTD_idx structure containing 
%           >> elec: {nx1 cell}
%           >> location: {nx1 cell}
%           >> offset: 2
%           >> nb_Gpix [nx1 double]
%           >> nb_Wpix [nx1 double]
%           >> PTD_idx.PTD [nx1 double]
%
% with:
% - elec: electrodes names
% - location: brain region where the centroid of each electrode is localized based on freesurfer parcellation
% - Offset: max # of voxels from contact centroid to include in PTD computation. 2 produces a 3x3 cube around
%   each contact's center voxel. correspond to the size of the cube around each electrode used to approximate 
%   the PTD (default = 2)
% - nb_Gpix and nb_Wpix correspond to the number of neocortical gray or white matter voxels within a cube centered around the electrode centroid
%   (default size = 3x3 ;centroid plus offset)
% - PTD proximal tissue density of white and neocortical gray matter around the electrode
%
% In addition the following files are created in the elec_recon subfolder
% of the patient's FreeSurfer directory:
%   -GreyWhite_classifications.mat
%   -GreyWhite_classifications.txt
%
% Function dependency:
% - MRIread from freesurfer (https://surfer.nmr.mgh.harvard.edu/)
%
% Files needed:
% - MRI from elec_recon (aparc+aseg.mgz)
% - electrodes coordinates (*.LEPTOVOX)
% - electrodes names (*.electrodeNames)
% - parcellation code table (FreeSurferColorLUT.txt)
%
%
% Reference:
% Mercier, M. R., Bickel, S., Megevand, P., Groppe, D. M., Schroeder, C. E.,
% Mehta, A. D., & Lado, F. A. (2017). Evaluation of cortical local field
% potential diffusion in stereotactic electro-encephalography recordings:
% A glimpse on white matter signal. NeuroImage, 147, 219-232.

% Change Log:
% 08-2017: a few other minor changes for iELVis by DG. In particular,
% instead of using wmparc.mgz, we now use aparc+aseg.mgz.
% 08-2017: adapted for iElvis by MrM;
% 02-2016: created by MrM;
%

% load parcellation file
fs_dir=getFsurfSubDir();
recon_folder=fullfile(fs_dir,fs_subj,'elec_recon');
parc_file=fullfile(fs_dir,fs_subj,'mri','aparc+aseg.mgz');
mri=MRIread(parc_file);

% load electrodes name
files=dir(fullfile(recon_folder,'*.electrodeNames'));
n=size(files,1);
if n==1
    label_file=fullfile(recon_folder,files.name);
elseif n==0
    disp('No electrodeNames file found. Please do it manualy.');
    [temp_file,elec_dir]=uigetfile(fullfile(recon_folder,'*.electrodeNames'),'select electrode names file');
    label_file=fullfile(elec_dir,temp_file);
    clear elec_dir temp_file files n
elseif n>1
    disp('More than one electrodeNames file found. Please do it manualy.');
    [temp_file,elec_dir]=uigetfile(fullfile(recon_folder,'*.electrodeNames'),'select electrode names file');
    label_file=fullfile(elec_dir,temp_file);
    clear elec_dir temp_file files n
end
elecInfo=csv2Cell(label_file,' ',2);
nElec=size(elecInfo,1);
label(1:nElec,1)={''};
for i=1:nElec
    label{i,1}=strcat(elecInfo{i,1},'_',elecInfo{i,2},'_',elecInfo{i,3});
end
clear elecInfo

% load electrodes coordinate
files=dir(fullfile(recon_folder,'*.LEPTOVOX'));
n=size(files,1);
if n==1
    elec_file=fullfile(recon_folder,files.name);
elseif n==0
    disp('No *.LEPTOVOX file found. Please do it manualy');
    [temp_file,elec_dir]=uigetfile(fullfile(recon_folder,'*.LEPTOVOX'),'select electrode coordinate file');
    elec_file=fullfile(elec_dir,temp_file);
    clear elec_dir temp_file files n
elseif n>1
    disp('More than one *.LEPTOVOX file found. Please choose it manualy');
    [temp_file,elec_dir]=uigetfile(fullfile(recon_folder,'*.LEPTOVOX'),'select electrode coordinate file');
    elec_file=fullfile(elec_dir,temp_file);
    clear elec_dir temp_file files n
end
coordCsv=csv2Cell(elec_file,' ',2);
elec=zeros(nElec,3);
for a=1:nElec
    for b=1:3
        elec(a,b)=str2double(coordCsv{a,b});
    end
end
clear coordCsv
elec=elec+1; % Voxel indexing starts at 0 but MATLAB indexing starts with 1

%% load look-up table for the FreeSurfer MRI atlases
FS_color_file = which('FreeSurferColorLUTnoFormat.txt');
if isempty(FS_color_file)
    disp('The freesurfer color code file, FreeSurferColorLUTnoFormat.txt, was not found. Please select it manualy.');
    [temp_file,elec_dir]=uigetfile(fullfile(recon_folder,'*.txt'),'Select freesurfer color code file');
    FS_color_file=fullfile(elec_dir,temp_file);
    clear elec_dir temp_file
else
    fprintf('Loading file %s\n',FS_color_file);
end

fid=fopen(FS_color_file);
C=textscan(fid,'%d %s %d %d %d %d');
fclose(fid);

region_lookup=C{:,2};
region_codes=C{:,1};
clear fid C p

%% find the proportion of neocortical grey and white matter surrounding the electrodes
sVol=size(mri.vol);
offset = 2;
for e=1:size(elec,1)
    disp(['Finding amount of surrounding grey and white matter for channel ' label{e}]);
    x=round(elec(e,2)); % switches x and y axes to match FreeSurfer
    y=round(elec(e,1)); % switches x and y axes to match FreeSurfer
    z=round(sVol(3)-elec(e,3)); %need to flip direction of last coordinate
    
    % DG Code for double checking electrode coordinates
    if e==0, %68, 61
        mri_mn=min(min(min(mri.vol)));
        mri_mx=max(max(max(mri.vol)));
        
        figure(2)
        colormap('gray');
        clf;
        subplot(1,3,1);
        imagesc(squeeze(mri.vol(:,y,:)),[mri_mn, mri_mx]);
        hold on;
        axis square;
        set(gca,'xdir','reverse');
        plot(z,x,'g.');
        
        subplot(1,3,2);
        imagesc(squeeze(mri.vol(x,:,:)),[mri_mn, mri_mx]);
        hold on;
        axis square;
        set(gca,'xdir','reverse');
        plot(z,y,'g.');
        
        subplot(1,3,3);
        imagesc(squeeze(mri.vol(:,:,z)),[mri_mn, mri_mx]);
        hold on;
        axis square;
        plot(y,x,'g.');
        
        disp('Done.');
    end
    
    % original labelling from parcellation
    ROI(e,1) =region_lookup(region_codes==mri.vol(x,y,z));
    
    % check that the contact is in gray or white matter
    gray_white=[strfind(lower(ROI{e,1}),'ctx') strfind(lower(ROI{e,1}),'cortex') ...
        strfind(lower(ROI{e,1}),'wm') strfind(lower(ROI{e,1}),'white')];
    if ~isempty(gray_white),
        % get the euclidean distances between the electrode and every voxel
        % in the MRI (this could be speed up a lot)
        for i=1:mri.volsize(1)
            for j=1:mri.volsize(2)
                for k=1:mri.volsize(3)
                    distances(i,j,k)=sqrt((i-x)^2+(j-y)^2+(k-z)^2); % could be replaced by ~ norm([i j k] - [x y z])
                end
            end
        end
        
        % do not include the offset
        tmp=mri.vol(distances<offset);
        % find the regions
        tmp_region ={};
        for i = 1:length(tmp)
            tmp_region(i,1) =region_lookup(region_codes==tmp(i));
        end
        % find gray matter voxel in the vicinity
        gm_nb = length(cell2mat(strfind(lower(tmp_region),'ctx')));
        gm_nb = gm_nb + length(cell2mat(strfind(lower(tmp_region),'cortex')));
        gm_w = gm_nb;
        % find white matter voxel in the vicinity
        wm_nb = length(cell2mat(strfind(lower(tmp_region),'wm')));
        wm_nb = wm_nb + length(cell2mat(strfind(lower(tmp_region),'white')));
        wm_w = wm_nb;
        
        ROI{e,2} = gm_w;
        ROI{e,3} = wm_w;
    else
        warning(['channel ' label{e} ' does not have its centroid in Neocortical Gray or White matter >> PTD=NaN']);
        ROI{e,2} = NaN;
        ROI{e,3} = NaN;
    end
end
%% write output file

for i=1:length(ROI)
    PTD_idx.elec(i,1)     = label(i);
    PTD_idx.location(i,1) = ROI(i,1);
    PTD_idx.nb_Gpix(i,1)  = cell2mat(ROI(i,2));
    PTD_idx.nb_Wpix(i,1)  = cell2mat(ROI(i,3));
    PTD_idx.PTD (i,1)     = (cell2mat(ROI(i,2)) - cell2mat(ROI(i,3))) / (cell2mat(ROI(i,2)) + cell2mat(ROI(i,3)));
    if (cell2mat(ROI(i,2)) + cell2mat(ROI(i,3))) ~= power(offset+1,3)
        warning(['channel ' label{e} ' has in its surrounding voxels that are neither labelled Gray or White matter; ' char(10)...
            'those voxels were not taking into account in PTD computation (see field nb_Gpix and nb_Wpix in the output)']);
    end
    % % otherwise a less strict version of the PTD taking into account the surrounding voxels that do not belong to Gray or White matter
    %     PTD_idx.PTD (i,1)     = (cell2mat(ROI(i,2)) - cell2mat(ROI(i,3))) / power(offset+1,3);
end
PTD_idx.offset = offset;

% Save PTD info as mat file
outfile_mat=fullfile(recon_folder,'GreyWhite_classifications.mat');
fprintf('Writing PTD values to %s\n',outfile_mat);
save(outfile_mat,'PTD_idx');

% Raw output in txt
outfile=fullfile(recon_folder,'GreyWhite_classifications.txt');
fprintf('Writing PTD values to %s\n',outfile);
fid=fopen(outfile,'w');
fprintf(fid,'%s\t %s\t %s\t %s\r\n','Electrode','location','% Grey with offset=2','% White with offset=2');
for i=1:size(ROI,1)
    fprintf(fid,'%s\t %s\t %d\t %d\t \r\n',label{i},ROI{i,:});
end
fclose(fid);
end