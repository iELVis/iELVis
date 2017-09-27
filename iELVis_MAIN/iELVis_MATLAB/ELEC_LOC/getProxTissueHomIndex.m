function proxTissueInfo = getProxTissueHomIndex(fsSub)
% function proxTissueInfo = getProxTissueHomIndex(fsSub)
%
% Computes Proxim Tissue Homogeneity (PTD) for each electrode.
%
% PTD is an index reflecting how anatomically homogenous the region around
% each electrode is in its local neighborhood (i.e., the 3x3 voxel cube
% centered on the contact). It varies from 1 (all voxels in the 3x3 cube
% are the same type of brain region) to 0.037 (the center voxel in the 3x3
% cube differs from all the others).
%
%
% Input: 
%   fsSub - name of the patient's FreeSurfer subject folder
%
% Output: 
%   proxTissueInfo structure containing 
%           elec_label: electrode name
%           pptn_homog: electrode's PTD value
%           anat_label: label of the anatomical region at the electrode's
%                       center
%
%
% Author: David Groppe
% Krembil Neuroscience Center (Sept. 2017)

% load parcellation file
fs_dir=getFsurfSubDir();
recon_folder=fullfile(fs_dir,fsSub,'elec_recon');
parc_file=fullfile(fs_dir,fsSub,'mri','aparc+aseg.mgz');
mri=MRIread(parc_file);

% load electrodes name
files=dir(fullfile(recon_folder,'*.electrodeNames'));
n=size(files,1);
if n==1
    label_file=fullfile(recon_folder,files.name);
elseif n==0
    disp('No electrodeNames file found. Please do it manually.');
    [temp_file,elec_dir]=uigetfile(fullfile(recon_folder,'*.electrodeNames'),'select electrode names file');
    label_file=fullfile(elec_dir,temp_file);
    clear elec_dir temp_file files n
elseif n>1
    disp('More than one electrodeNames file found. Please do it manually.');
    [temp_file,elec_dir]=uigetfile(fullfile(recon_folder,'*.electrodeNames'),'select electrode names file');
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
files=dir(fullfile(recon_folder,'*.LEPTOVOX'));
n=size(files,1);
if n==1
    elec_file=fullfile(recon_folder,files.name);
elseif n==0
    disp('No *.LEPTOVOX file found. Please do it manually');
    [temp_file,elec_dir]=uigetfile(fullfile(recon_folder,'*.LEPTOVOX'),'select electrode coordinate file');
    elec_file=fullfile(elec_dir,temp_file);
    clear elec_dir temp_file files n
elseif n>1
    disp('More than one *.LEPTOVOX file found. Please choose it manually');
    [temp_file,elec_dir]=uigetfile(fullfile(recon_folder,'*.LEPTOVOX'),'select electrode coordinate file');
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
elec=elec+1; % Voxel indexing starts at 0 but MATLAB indexing starts with 1

%% load look-up table for the FreeSurfer MRI atlases
% FS_color_file = which('FreeSurferColorLUTnoFormat.txt');
% if isempty(FS_color_file)
%     disp('The freesurfer color code file, FreeSurferColorLUTnoFormat.txt, was not found. Please select it manually.');
%     [temp_file,elec_dir]=uigetfile(fullfile(recon_folder,'*.txt'),'Select freesurfer color code file');
%     FS_color_file=fullfile(elec_dir,temp_file);
%     clear elec_dir temp_file
% else
%     fprintf('Loading file %s\n',FS_color_file);
% end
% 
% fid=fopen(FS_color_file);
% C=textscan(fid,'%d %s %d %d %d %d');
% fclose(fid);
% 
% region_lookup=C{:,2};
% region_codes=C{:,1};
% clear fid C p

% DG
% asegFname=fullfile(fs_dir,fsSub,'mri','aparc+aseg.mgz');
% %asegFname=[fs_dir '/' fsSub '/mri/aparc.a2009s+aseg.mgz'];
% if ~exist(asegFname,'file')
%    error('File %s not found.',asegFname); 
% end
% 
% aseg=MRIread(asegFname);

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
% id=find(aseg.vol(coordILA(1),coordILA(2),coordILA(3))==tbl{1});
% id=min(id);
% anatLabel=tbl{2}{id};

%% find the proportion of neocortical grey and white matter surrounding the electrodes
sVol=size(mri.vol);
offset = 2;
n_elec=size(elec,1);
pptn_homogeneous=zeros(n_elec,1);
anatLabel=cell(n_elec,1);
for e=1:n_elec,
    disp(['Finding amount of surrounding grey and white matter for channel ' label{e}]);
    x=round(elec(e,2)); % switches x and y axes to match FreeSurfer
    y=round(elec(e,1)); % switches x and y axes to match FreeSurfer
    z=round(sVol(3)-elec(e,3)); %need to flip direction of last coordinate
    
    id=find(mri.vol(x,y,z)==tbl{1});
    id=min(id);
    anatLabel{e}=tbl{2}{id};

    
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
    %ROI(e,1) =region_lookup(region_codes==mri.vol(x,y,z));

    % get the euclidean distances between the electrode and every voxel
    % in the MRI (this could be sped up a lot)
    nbor_ct=0;
    for i=1:mri.volsize(1)
        for j=1:mri.volsize(2)
            for k=1:mri.volsize(3)
                temp_dst=sqrt((i-x)^2+(j-y)^2+(k-z)^2); % could be replaced by ~ norm([i j k] - [x y z])
                if temp_dst<offset,
                   nbor_ct=nbor_ct+1;
                   if mri.vol(x,y,z)==mri.vol(i,j,k),
                       pptn_homogeneous(e)=pptn_homogeneous(e)+1;
                   end
                end
            end
        end
    end
    pptn_homogeneous(e)=pptn_homogeneous(e)/nbor_ct;
end

%% Clean up labels
elec_labels=cell(n_elec,1);
for eloop=1:n_elec,
    split_id=find(label{eloop}=='_');
    elec_labels{eloop}=label{eloop}(1:split_id-1);
end

%% Add rest of info to struct
proxTissueInfo=struct('elec_label',elec_labels);
for eloop=1:n_elec,
    proxTissueInfo(eloop).pptn_homog=pptn_homogeneous(eloop);
    proxTissueInfo(eloop).anat_label=anatLabel{eloop};
end

%% write output file

% for i=1:length(ROI)
%     PTD_idx.elec(i,1)     = label(i);
%     PTD_idx.location(i,1) = ROI(i,1);
%     PTD_idx.nb_Gpix(i,1)  = cell2mat(ROI(i,2));
%     PTD_idx.nb_Wpix(i,1)  = cell2mat(ROI(i,3));
%     PTD_idx.PTD (i,1)     = (cell2mat(ROI(i,2)) - cell2mat(ROI(i,3))) / (cell2mat(ROI(i,2)) + cell2mat(ROI(i,3)));
%     if (cell2mat(ROI(i,2)) + cell2mat(ROI(i,3))) ~= power(offset+1,3)
%         warning(['channel ' label{e} ' has in its surrounding voxels that are neither labelled Gray or White matter; ' char(10)...
%             'those voxels were not taking into account in PTD computation (see field nb_Gpix and nb_Wpix in the output)']);
%     end
%     % % otherwise a less strict version of the PTD taking into account the surrounding voxels that do not belong to Gray or White matter
%     %     PTD_idx.PTD (i,1)     = (cell2mat(ROI(i,2)) - cell2mat(ROI(i,3))) / power(offset+1,3);
% end
% PTD_idx.offset = offset;

% Save PTD info as mat file
% outfile_mat=fullfile(recon_folder,'GreyWhite_classifications.mat');
% fprintf('Writing PTD values to %s\n',outfile_mat);
% save(outfile_mat,'PTD_idx');
% 
% % Raw output in txt
% outfile=fullfile(recon_folder,'GreyWhite_classifications.txt');
% fprintf('Writing PTD values to %s\n',outfile);
% fid=fopen(outfile,'w');
% fprintf(fid,'%s\t %s\t %s\t %s\r\n','Electrode','location','% Grey with offset=2','% White with offset=2');
% for i=1:size(ROI,1)
%     fprintf(fid,'%s\t %s\t %d\t %d\t \r\n',label{i},ROI{i,:});
% end
% fclose(fid);
end