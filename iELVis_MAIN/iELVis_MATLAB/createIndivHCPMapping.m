function createIndivHCP(fsSub)
% function createIndivYeoMapping(fsSub)
%
% This function assigns each point on an individual subject's pial surface
% to the HCP atlas. It simply takes the mapping of the individual brain to
% that of the FreeSurfer average brain and assigns each point in the
% individual the label of the closest point in the average brain.
%
% Input:
%   fsSub = FreeSurfer fsSubect name
%
% Outputs:
%   The following MATLAB files are created in the subject's label
%   FreeSurfer subdirectory:
%    lh.HCP-MMP1.annot
%    rh.HCP-MMP1.annot
%
% Author:
% Noah Markowitz
% The Human Brain Mapping Labratory
% Feinstein Institutes for Medical Research
% October 2021
%

% Get location of FreeSurfer directories
fsDir=getFsurfSubDir();
labelFolder=fullfile(fsDir,fsSub,'label');

avg_dir=fullfile(fsDir,'fsaverage');
sub_dir=fullfile(fsDir,fsSub);

% Check if the HCP atlas files have been downloaded
if ~exist( fullfile(avg_dir, 'label', 'lh.HCP-MMP1.annot'), 'file')
    downloadHCPatlas(1);
end

for hemLoop=1:2,
    if hemLoop==1,
        hem='lh';
    else
        hem='rh';
    end
    fprintf('Creating HCP atlas for hemisphere: %s\n',hem);
    
    
    %% Read Sub Pial Surf
    fname=[sub_dir filesep 'surf' filesep hem '.pial'];
    pial=read_surf_helper(fname);
    
    %% Read Sub Spherical Surf
    fname=[sub_dir filesep 'surf' filesep hem '.sphere.reg'];
    sph=read_surf_helper(fname);
    n_sub_vert=size(sph.vert,1);
    
    
    %% Load Avg Spherical Surf
    fname=[avg_dir filesep 'surf' filesep hem '.sphere.reg'];
    avg_sph=read_surf_helper(fname);
    n_avg_vert=length(avg_sph.vert);
    
    
    %% Load HCP atlas
    fnameHCP = [hem '.HCP-MMP1.annot'];
    [avgBrain, label, colortable]=read_annotation(fullfile(fsDir,'fsaverage','label',fnameHCP));
    
    indivBrainHCP=zeros(n_sub_vert,1);
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
        indivBrainHCP(b)=label(id);
    end
    
    %% Write out the annotation file
    
    annotFname=fullfile(labelFolder,[hem '.HCP-MMP1.annot']);
    fprintf('Saving HCP-MMP1 %s\n',annotFname);    
    write_annotation(annotFname, [0:(n_sub_vert-1)] ,indivBrainHCP, colortable);
    
end
