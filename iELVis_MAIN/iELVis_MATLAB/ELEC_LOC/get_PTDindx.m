function WandG = get_PTDindx(fs_subj)
%
% compute Proximal Tissue Density for each electrode
% PTD is an index reflect the density of gray and white matter surrounding
% a stereotactic electrode implanted in the brain.
% as describe in: Mercier et al., Neuroimage 2017
% 
% Please note that this function belong to the iElvis toolbox
% and is therefore subjected to the same regulations
% (contacts: david.m.groppe@gmail.com ; manuel.mercier@a3.epfl.ch)
%
% input: fs_subj name of the subject in the FS folder
%
% output: WandG strcture containing 
%           >> elec: {nx1 cell} (from 
%           >> location: {nx1 cell}
%           >> offset: 2
%           >> nb_Gpix [nx1 double]
%           >> nb_Wpix [nx1 double]
%           >> WandG.PTD [nx1 double]
%
% with:
% - Offset correspond to the size of the cube (default = 2)
% - nb_Gpix and nb_Wpix correspond to the number of Gray or White matter voxels
%   within a cube centered around the electrode centroid
%   (default size = 27 ;centroid plus offset)
% - PTD proximal tissue density of gray and White matter electrode around the electrode
%
%
% functions needed:
% - MRIread from freesurfer
% - euclid_dist_mtx from epiproc
% 
% files needed:
% - MRI from elec_recon (wmparc.mgz)
% - electrodes coordinates (*.3dUndump.VOX)
% - electrodes names (*.electrodeNames)
% - parcellation code table (FreeSurferColorLUT.txt)
% - information about channel (*_channelinfo.mat)
% 
% if hemispheres were processed separatlty,
% use the last section of the script to merge
% files from Left and Right infos ...
%
%
% 08-2017: addapted for iElvis by MrM;
% 02-2016: created by MrM;
%


%load parcellation file
if isunix || ismac
fs_dir=getenv('SUBJECTS_DIR');
recon_folder=fullfile(fs_dir,fs_subj,'elec_recon/');
parc_file=fullfile(fs_dir,fs_subj,'mri','wmparc.mgz');
mri=MRIread(parc_file);
elseif ispc
fs_dir=uigetdir(userpath,'select freesurfer directory');
recon_folder=fullfile(fs_dir,fs_subj,'elec_recon\');
parc_file=fullfile(fs_dir,fs_subj,'mri','wmparc.mgz');
mri=MRIread([parc_file(1:end-4) '.mgz']);
end        
%load coordinates
% files=dir([recon_folder '*.3dUndump.VOX']);
files=dir([recon_folder '*.LEPTOVOX']);

n=size(files,1);
if n==1
    elec_file=fullfile(recon_folder,files.name);
elseif n==0
    disp('No 3dUndump file found. Please do it manualy');
    [temp_file,elec_dir]=uigetfile([recon_folder '*.LEPTOVOX'],'select electrode coordinate file');
    elec_file=fullfile(elec_dir,temp_file);    
    clear elec_dir temp_file    
elseif n>1
    disp('More than one 3dUndumpfile found. Please choose it manualy');
    [temp_file,elec_dir]=uigetfile([recon_folder '*.LEPTOVOX'],'select electrode coordinate file');
    elec_file=fullfile(elec_dir,temp_file);
    clear elec_dir temp_file
end
%load electrodes names
files=dir([recon_folder '*.electrodeNames']);
n=size(files,1);
if n==1
    label_file=fullfile(recon_folder,files.name);
elseif n==0
    disp('No electrodeNames file found. Please do it manualy.');
    [temp_file,elec_dir]=uigetfile([recon_folder '*.electrodeNames'],'select electrode names file');
    label_file=fullfile(elec_dir,temp_file);
    clear elec_dir temp_file    
elseif n>1
    disp('More than one electrodeNames file found. Please do it manualy.');
    [temp_file,elec_dir]=uigetfile([recon_folder '*.electrodeNames'],'select electrode names file');
    label_file=fullfile(elec_dir,temp_file);
    clear elec_dir temp_file
end

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
%% load electrode coordinates
fid=fopen(elec_file);
tmp=textscan(fid,'%s %s %s');
fclose(fid);
for i=3:length(tmp{1})
    elec(i-2,1)=str2num(tmp{1}{i});
    elec(i-2,2)=str2num(tmp{2}{i});
    elec(i-2,3)=str2num(tmp{3}{i});
end
clear tmp;

% load electrode labels
fid=fopen(label_file);
tmp=textscan(fid,'%s %s %s');
fclose(fid);
for i=3:length(tmp{1})
    label{i-2,1}=strcat(tmp{1}{i},'_',tmp{2}{i},'_',tmp{3}{i});
end
clear tmp;



%% find the proportion of grey and white matter surrounding the electrodes

offset = 2;
for e=1:size(elec,1)
    disp(['Finding amount of surrounding grey and white matter for channel ' label{e}]);
    x=round(elec(e,2)); % switches x and y planes to match FreeSurfer
    y=round(elec(e,1)); % switches x and y planes to match FreeSurfer
    z=round(elec(e,3));
    
    % original labelling from parcellation
    ROI(e,1) =region_lookup(region_codes==mri.vol(x,y,z));
    % get the euclidean distances between the electrode and every voxel in the MRI
    for i=1:mri.volsize(1)
        for j=1:mri.volsize(2)
            for k=1:mri.volsize(3)
                distances(i,j,k)=sqrt((i-x)^2+(j-y)^2+(k-z)^2);
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
    
    ROI{e,2*(offset-1)} = gm_w;
    ROI{e,3*(offset-1)} = wm_w;
end
%% write output file

for i=1:length(ROI)   
    WandG.elec(i,1)     = label(i);
    WandG.location(i,1) = ROI(i,1);
    WandG.nb_Gpix(i,1)  = cell2mat(ROI(i,2));
    WandG.nb_Wpix(i,1)  = cell2mat(ROI(i,3));
    WandG.PTD (i,1)     = (cell2mat(ROI(i,2)) - cell2mat(ROI(i,3))) / (cell2mat(ROI(i,2)) + cell2mat(ROI(i,3)));
end
WandG.offset = offset;
%save
save([recon_folder 'GreyWhite_classifications.mat'],'WandG');
% raw output in txt
outfile=([recon_folder 'GreyWhite_classifications.txt']);
fid=fopen(outfile,'w');
fprintf(fid,'%s\t %s\t %s\t %s\r\n','Electrode','location','% Grey with offset=2','% White with offset=2');
for i=1:size(ROI,1)
    fprintf(fid,'%s\t %s\t %d\t %d\t \r\n',label{i},ROI{i,:});
end
fclose(fid);



%% if separate R/L files you can run this piece of code to merge the L/R outputfiles
% 
% for i = 1:length(WandG_L.elec)
%     if ~isempty(WandG_L.elec{i});
%         WandG.elec{i,1}       = WandG_L.elec{i};
%         WandG.location{i,1}   = WandG_L.location{i};
%         WandG.G1(i,1)       = WandG_L.G1(i);
%         WandG.W1(i,1)       = WandG_L.W1(i);
% %         WandG.G2(i,1)       = WandG_L.G2(i);
% %         WandG.W2(i,1)       = WandG_L.W2(i);
% %         WandG.G3(i,1)       = WandG_L.G3(i);
% %         WandG.W3(i,1)       = WandG_L.W3(i); 
%     elseif isempty(WandG_L.elec{i}) && i <= length(WandG_R.elec);
%         WandG.elec{i,1}       = WandG_R.elec{i};
%         WandG.location{i,1}   = WandG_R.location{i};
%         WandG.G1(i,1)       = WandG_R.G1(i);
%         WandG.W1(i,1)       = WandG_R.W1(i);
% %         WandG.G2(i,1)       = WandG_R.G2(i);
% %         WandG.W2(i,1)       = WandG_R.W2(i);
% %         WandG.G3(i,1)       = WandG_R.G3(i);
% %         WandG.W3(i,1)       = WandG_R.W3(i); 
%     end
% end
end