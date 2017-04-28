% function ntools_elec_CTautocoreg
clear all; close all
% using flirt dof 6 to coregister the postop CT to preop MRI (usually FS
% recon T1)

%coreg = '/home/halgdev/projects/nyuproj/loc/'; % your localization folder
coreg= '/Users/dgroppe/LOC';
subj = getenv('SUBJECTS_DIR');

t1 = menu('Select the pre-operation image','Pick up my own T1','Select Freesurfer T1.mgz');
if t1==1
    [preop, preoppath] = uigetfile('*.*', 'Select the pre-operation image',coreg);
else
    [preop, preoppath] = uigetfile('*.*', 'Select the pre-operation image',subj);
end
preop = fullfile(preoppath,preop);

[elec, elecpath] = uigetfile('*.*', 'Select the post-operation image with electrodes',coreg);
elec = fullfile(elecpath,elec);

% convert format and orientation
% preop
[~,name_preop,ext] = fileparts(preop);
[~,name_preop] = fileparts(name_preop);
preop_nii = [elecpath,name_preop,ext];
preop_nii = regexprep(preop_nii,ext,'.nii.gz');
cmd_convert = sprintf('mri_convert --out_orientation RAS %s %s',preop,preop_nii);
[status, msg] = unix(cmd_convert);
if status, disp(msg); return; end

% elec
elec_nii = regexprep(elec,'.img','.nii.gz');
cmd_convert = sprintf('mri_convert --out_orientation RAS %s %s',elec,elec_nii);
[status, msg] = unix(cmd_convert);
if status, disp(msg); return; end

% coregister
[~,name_elec] = fileparts(elec_nii);
[~,name_elec] = fileparts(name_elec);
elec_preop = [elecpath,name_elec,'_preop.nii.gz'];
elec_preop_brain =  [elecpath,name_elec,'_preop_brain.nii.gz'];
elec_preop_omat = [elecpath,name_elec,'_preop.mat'];
unix(sprintf('flirt -in %s -ref %s -out %s -omat %s -dof 6 -cost normmi -interp nearestneighbour',elec_nii,preop_nii,elec_preop,elec_preop_omat),'-echo');
unix(sprintf('bet %s %s_brain.nii.gz -f 0.5 -g 0 -m',preop_nii,[elecpath,name_preop]),'-echo');
unix(sprintf('fslmaths %s -mas %s_brain_mask.nii.gz %s',elec_preop,[elecpath,name_preop],elec_preop_brain));

% check results
unix(sprintf('fslview %s %s',preop_nii,elec_preop));


