function createIndivYeoMapping(fsSub)
% function createIndivYeoMapping(fsSub)
%
% This function assigns each point on an individual subject's pial surface
% to the Yeo-7 area and Yeo-17 area atlases, which are based on resting
% state fMRI data. It simply takes the mapping of the individual brain to
% that of the FreeSurfer average brain and assigns each point in the
% individual the label of the closest point in the average brain.
%
% Input:
%   fsSub = FreeSurfer fsSubect name
%
% Outputs:
%   The following MATLAB files are created in the subject's label
%   FreeSurfer subdirectory:
%    lh_Yeo2011_17Networks_N1000.mat
%    rh_Yeo2011_17Networks_N1000.mat
%    lh_Yeo2011_7Networks_N1000.mat
%    rh_Yeo2011_7Networks_N1000.mat
%
% Author:
% David Groppe
% Honeylab
% January, 2016
%

% Get location of FreeSurfer directories
fsDir=getFsurfSubDir();
labelFolder=fullfile(fsDir,fsSub,'label');

avg_dir=[fsDir '/' 'fsaverage'];
sub_dir=[fsDir '/' fsSub];

%% 7 area labels taken from the original paper:
% Yeo BT, Krienen FM, Sepulcre J, Sabuncu MR, Lashkari D, Hollinshead M, 
% Roffman JL, Smoller JW, Zollei L., Polimeni JR, Fischl B, Liu H, Buckner 
% RL. The organization of the human cerebral cortex estimated by intrinsic 
% functional connectivity. J Neurophysiol 106(3):1125-65, 2011.

y7labels=cell(1,7);
y7labels{1}='MedialWall';
y7labels{2}='Visual';
y7labels{3}='Somatomotor';
y7labels{4}='Dorsal Attention';
y7labels{5}='Ventral Attention';
y7labels{6}='Limbic';
y7labels{7}='Frontoparietal';
y7labels{8}='Default';

for hemLoop=1:2,
    if hemLoop==1,
        hem='lh';
    else
        hem='rh';
    end
    fprintf('Creating Yeo mapping for hemisphere: %s\n',hem);
    
    
    %% Read Sub Pial Surf
    fname=[sub_dir '/surf/' hem '.pial'];
    pial=read_surf_helper(fname);
    
    
    %% Read Sub Spherical Surf
    fname=[sub_dir '/surf/' hem '.sphere.reg'];
    sph=read_surf_helper(fname);
    n_sub_vert=size(sph.vert,1);
    
    
    %% Load Avg Spherical Surf
    fname=[avg_dir '/surf/' hem '.sphere.reg'];
    avg_sph=read_surf_helper(fname);
    n_avg_vert=length(avg_sph.vert);
    
    
    %% Load Yeo atlases
    fname7=[hem '.Yeo2011_7Networks_N1000.annot'];
    fname17=[hem '.Yeo2011_17Networks_N1000.annot'];
    [avgBrainYeo7, label7, colortable7]=read_annotation(fullfile(fsDir,'fsaverage','label',fname7));
    [avgBrainYeo17, label17, colortable17]=read_annotation(fullfile(fsDir,'fsaverage','label',fname17));
    
    for b=2:8,
        colortable7.struct_names{b}=y7labels{b};
    end
    
    indivBrainYeo7=zeros(n_sub_vert,1);
    indivBrainYeo17=zeros(n_sub_vert,1);
    vertices=zeros(n_sub_vert,1);
    
    
    %% Map pial surface vertices in subject's sphere to avg sph
    fprintf('Processing vertex:\n');
    for b=1:n_sub_vert,
        if ~rem(b,1000)
            fprintf('%d of %d\n',b,n_sub_vert);
        end
        dst=sum( (avg_sph.vert-repmat(sph.vert(b,:),n_avg_vert,1)).^2 ,2);
        [dummy id]=min(dst);
        
        vertices(b)=id;
        indivBrainYeo7(b)=label7(id);
        indivBrainYeo17(b)=label17(id);
    end
    
    %% Export individual annotation as a mat file
    annotFname7=fullfile(labelFolder,[hem '_Yeo2011_7Networks_N1000.mat']);
    fprintf('Saving Yeo7 %s\n',annotFname7);
    label=indivBrainYeo7;
    colortable=colortable7;
    save(annotFname7,'label','colortable');
    
    annotFname17=fullfile(labelFolder,[hem '_Yeo2011_17Networks_N1000.mat']);
    fprintf('Saving Yeo17 %s\n',annotFname17);
    label=indivBrainYeo17;
    colortable=colortable17;
    save(annotFname17,'label','colortable');
    
end
