% function ntools_elec_autocoreg
clear all; close all
% using flirt dof 6 to coregister the postop MRI to preop MRI (usually FS
% recon T1)

coreg = pwd; % your localization folder
subj = getenv('SUBJECTS_DIR');

t1 = menu('Select the pre-operation image','Pick up my own T1','Select Freesurfer T1.mgz');
if t1==1
    [preop, preoppath] = uigetfile('*.*', 'Select the pre-operation image',coreg);
else
    [preop, preoppath] = uigetfile('*.*', 'Select the pre-operation image',subj);
end
preop = fullfile(preoppath,preop);

[aseg,asegpath] = uigetfile('*.mgz','Select subject aseg file',preoppath);
aseg = fullfile(asegpath,aseg);

[elec, elecpath] = uigetfile('*.*', 'Select the post-operation image with electrodes',coreg);
elec = fullfile(elecpath,elec);

%% convert format and orientation
% preop
[~,name_preop,ext] = fileparts(preop);
[~,name_preop] = fileparts(name_preop);
preop_nii = [elecpath,name_preop,ext];
preop_nii = regexprep(preop_nii,ext,'.nii.gz');
cmd_convert = sprintf('mri_convert --out_orientation RAS %s %s',preop,preop_nii);
[status, msg] = unix(cmd_convert);
if status, disp(msg); return; end

%% aseg
[~,name_aseg,ext] = fileparts(aseg);
[~,name_aseg] = fileparts(name_aseg);
aseg_nii = [elecpath,name_aseg,ext];
aseg_nii = regexprep(aseg_nii,ext,'.nii.gz');
cmd_convert = sprintf('mri_convert --out_orientation RAS %s %s',aseg,aseg_nii);
[status, msg] = unix(cmd_convert);
if status, disp(msg); return; end

%% elec
elec_nii = regexprep(elec,'.img','.nii.gz');
cmd_convert = sprintf('mri_convert --out_orientation RAS %s %s',elec,elec_nii);
[status, msg] = unix(cmd_convert);
if status, disp(msg); return; end

%% coregister
[~,name_elec] = fileparts(elec_nii);
[~,name_elec] = fileparts(name_elec);
elec_preop = [elecpath,name_elec,'_preop.nii.gz'];
elec_preop_brain =  [elecpath,name_elec,'_preop_brain.nii.gz'];
elec_preop_omat = [elecpath,name_elec,'_preop.mat'];
unix(sprintf('flirt -in %s -ref %s -out %s -omat %s -dof 6 -interp trilinear',elec_nii,preop_nii,elec_preop,elec_preop_omat),'-echo');
unix(sprintf('bet %s %s_brain.nii.gz -f 0.5 -g 0 -m',preop_nii,[elecpath,name_preop]),'-echo');
unix(sprintf('fslmaths %s -mul %s_brain_mask.nii.gz %s',elec_preop,[elecpath,name_preop],elec_preop_brain));

%% remove cerebellum
cerebellum = [elecpath,'cerebellum.nii.gz'];
unix(sprintf('mri_binarize --i %s --match 6 7 8 45 46 47 --o %s',aseg_nii,cerebellum));
% dialte cerebellum
unix(sprintf('fslmaths %s -dilM %s',cerebellum,cerebellum));
% inverse mask
cerebellum_inv = [elecpath,'cerebellum_inv.nii.gz'];
unix(sprintf('fslmaths %s -sub 1 -abs %s',cerebellum,cerebellum_inv));
elec_preop_cortex = regexprep(elec_preop_brain,'brain','cortex');
unix(sprintf('fslmaths %s -mas %s %s',elec_preop_brain,cerebellum_inv,elec_preop_cortex));

%% check results
unix(sprintf('fslview %s %s',elec_preop, elec_preop_cortex));

