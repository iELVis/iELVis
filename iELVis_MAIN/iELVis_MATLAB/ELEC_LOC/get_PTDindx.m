function PTD_idx = get_PTDindx(fs_subj)
%
% compute Proximal Tissue Density for each electrode
% as describe in Mercier et al., Neuroimage 2017
%
% PTD is an index reflecting the density of gray and white matter surrounding
% a stereotactic electrode that has it centroid either in the Gray or in the White matter.
% (otheriwse PTD = NaN)
% Be carefull when a contact is in the vicinity of subcortical structures
% PTD is computed exclusively by taking into account surrounding Gray and White matter voxels
% no other tissue will be taken into account (e.g. Hippocampus, Amygdala...)
% 
% Please note that this function belong to the iElvis toolboxe
% and is therefore subjected to the same regulations
% (contacts: david.m.groppe@gmail.com ; manuel.mercier@a3.epfl.ch)
%
% input: fs_subj name of the subject in the FS folder
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
% - elec: electrodes name
% - location: brain region where the centroid of each electrode is localized based on freesurfer parcellation
% - Offset: correspond to the size of the cube around each electrode used to approximate the PTD (default = 2)
% - nb_Gpix and nb_Wpix correspond to the number of Gray or White matter voxels within a cube centered around the electrode centroid
%   (default size = 27 ;centroid plus offset)
% - PTD proximal tissue density of gray and White matter electrode around the electrode
%
%
% functions dependency:
% - MRIread from freesurfer (https://surfer.nmr.mgh.harvard.edu/)
% 
% files needed:
% - MRI from elec_recon (wmparc.mgz)
% - electrodes coordinates (*.3dUndump.VOX)
% - electrodes names (*.electrodeNames)
% - parcellation code table (FreeSurferColorLUT.txt)
% 
%
%
% 08-2017: adapted for iElvis by MrM;
% 02-2016: created by MrM;
%


% load parcellation file
if isunix || ismac
fs_dir=getenv('SUBJECTS_DIR');
recon_folder=fullfile(fs_dir,fs_subj,'elec_recon/');
elseif ispc
fs_dir=uigetdir(userpath,'select freesurfer directory');
recon_folder=fullfile(fs_dir,fs_subj,'elec_recon\');
end
parc_file=fullfile(fs_dir,fs_subj,'mri','wmparc.mgz');
mri=MRIread(parc_file);

% load electrodes name
files=dir([recon_folder '*.electrodeNames']);
n=size(files,1);
if n==1
    label_file=fullfile(recon_folder,files.name);
elseif n==0
    disp('No electrodeNames file found. Please do it manualy.');
    [temp_file,elec_dir]=uigetfile([recon_folder '*.electrodeNames'],'select electrode names file');
    label_file=fullfile(elec_dir,temp_file);
    clear elec_dir temp_file files n
elseif n>1
    disp('More than one electrodeNames file found. Please do it manualy.');
    [temp_file,elec_dir]=uigetfile([recon_folder '*.electrodeNames'],'select electrode names file');
    label_file=fullfile(elec_dir,temp_file);
    clear elec_dir temp_file files n
end
fid=fopen(label_file);
tmp=textscan(fid,'%s %s %s');
fclose(fid);
for i=3:length(tmp{1})
    label{i-2,1}=strcat(tmp{1}{i},'_',tmp{2}{i},'_',tmp{3}{i});
end
clear tmp;

% load electrodes coordinate
% files=dir([recon_folder '*.3dUndump.VOX']);
files=dir([recon_folder '*.LEPTOVOX']);
n=size(files,1);
if n==1
    elec_file=fullfile(recon_folder,files.name);
elseif n==0
    disp('No 3dUndump file found. Please do it manualy');
    [temp_file,elec_dir]=uigetfile([recon_folder '*.LEPTOVOX'],'select electrode coordinate file');
    elec_file=fullfile(elec_dir,temp_file);    
    clear elec_dir temp_file files n
elseif n>1
    disp('More than one 3dUndumpfile found. Please choose it manualy');
    [temp_file,elec_dir]=uigetfile([recon_folder '*.LEPTOVOX'],'select electrode coordinate file');
    elec_file=fullfile(elec_dir,temp_file);
    clear elec_dir temp_file files n
end
fid=fopen(elec_file);
tmp=textscan(fid,'%s %s %s');
fclose(fid);
for i=3:length(tmp{1})
    elec(i-2,1)=str2num(tmp{1}{i});
    elec(i-2,2)=str2num(tmp{2}{i});
    elec(i-2,3)=str2num(tmp{3}{i});
end
clear tmp;

%% load look-up table for the FreeSurfer MRI atlases
FS_color_file = which('FreeSurferColorLUTnoformat.txt');
if isempty(FS_color_file)
    disp('the freesurfer color code file was not found. Please select it manualy.');
    [temp_file,elec_dir]=uigetfile([recon_folder '*.txt'],'select  freesurfer color code file');
    FS_color_file=fullfile(elec_dir,temp_file);
    clear elec_dir temp_file
end
    
fid=fopen(FS_color_file);
C=textscan(fid,'%d %s %d %d %d %d');
fclose(fid);

region_lookup=C{:,2};
region_codes=C{:,1};
clear fid C p

%% find the proportion of grey and white matter surrounding the electrodes

offset = 2;
for e=1:size(elec,1)
    disp(['Finding amount of surrounding grey and white matter for channel ' label{e}]);
    x=round(elec(e,2)); % switches x and y planes to match FreeSurfer
    y=round(elec(e,1)); % switches x and y planes to match FreeSurfer
    z=round(elec(e,3));
    
    % original labelling from parcellation
    ROI(e,1) =region_lookup(region_codes==mri.vol(x,y,z));
    
    % check that the contact is in gray or white matter
    if strfind(lower(ROI{e,1}),'ctx') || strfind(lower(ROI{e,1}),'cortex') || ...
       strfind(lower(ROI{e,1}),'wm') || strfind(lower(ROI{e,1}),'white')
    
        % get the euclidean distances between the electrode and every voxel in the MRI
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
        warning(['channel ' label{e} ' does not have its centroid in neither Gray or White matter >> PTD=NaN']);
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
%save
save([recon_folder 'GreyWhite_classifications.mat'],'PTD_idx');
% raw output in txt
outfile=([recon_folder 'GreyWhite_classifications.txt']);
fid=fopen(outfile,'w');
fprintf(fid,'%s\t %s\t %s\t %s\r\n','Electrode','location','% Grey with offset=2','% White with offset=2');
for i=1:size(ROI,1)
    fprintf(fid,'%s\t %s\t %d\t %d\t \r\n',label{i},ROI{i,:});
end
fclose(fid);
end