% Script for testing if iELVis is installed
%
%% Test to make sure iELVis in path
disp('Checking that iELVis files are in the MATLAB''s path...');
if ~exist('plotPialSurf','file')
    error('iELVis file plotPialSurf.m not found. You need to add iELVis folders to your MATLAB path.');
end


%% Test to make sure FreeSurfer subject directory exists
disp('Checking that FreeSurfer "subjects" directory is defined...');
try
    fsDir=getFsurfSubDir();
catch
    error('FreeSurfer "subjects" directory not found. Store the path in a variable called "globalFsDir" in your startup.m file.');
end


%% Test to make sure that PT001 is in path
disp('Checking that PT001 folder is in FreeSurfer "subjects" directory...');
if ~exist(fullfile(fsDir,'PT001'),'dir')
   error('Folder for PT001 is not present in FreeSurfer subjects directory. Download it from here https://osf.io/afwrz/ and add it.'); 
end


%% Test to make sure that fsaverage is in path
disp('Checking that fsaverage folder is in FreeSurfer "subjects" directory...');
if ~exist(fullfile(fsDir,'fsaverage'),'dir')
   error('Folder for fsaverage is not present in FreeSurfer subjects directory. Download it from here https://osf.io/qe7pz/ and add it.'); 
end



%% Test simplest brainshift correction
disp('Checking simplest brainshift correction method on demo data PT001...');
dykstraElecPjct('PT001',0);


%% Test electrode visualization with atlas
disp('Checking atlas visualization on demo data PT001...');
cfg=[];
cfg.view='l';
cfg.overlayParcellation='DK';
cfg.showLabels='y';
cfg.title='PT001: DK Atlas'; 
cfgOut=plotPialSurf('PT001',cfg);


%% Test electrode coloration
disp('Checking electrode coloration on demo data PT001...');
elecNames=cell(8,1);
for a=1:8,
    elecNames{a}=sprintf('LGd%d',a+16);
end
cfg=[];
cfg.view='l';
cfg.elecShape='sphere';
cfg.elecColors=rand(8,1);
cfg.elecColorScale='minmax';
cfg.showLabels='n';
cfg.elecUnits='r';
cfg.elecNames=elecNames;
cfg.elecSize=2;
cfg.title='PT001: Stimulus Correlations';
cfgOut=plotPialSurf('PT001',cfg);


%% Test neuroimaging overlay
disp('Checking fMRI visualization on demo data PT001...');
cfg=[];
cfg.view='lomni';
cfg.olayUnits='z';
cfg.pialOverlay=sprintf('%s/PT001/fmri/handMotorLH.mgh',fsDir);
cfgOut=plotPialSurf('PT001',cfg);


%% Test avg brain mapping
disp('Checking mapping to average brain on demo data PT001...');
[avgCoords, elecNames, isLeft]=sub2AvgBrain('PT001',[]);


%% Test assigning electrodes to atlas areas
disp('Checking mapping to electrodes to atlas areas using demo data PT001...');
parcOut=elec2Parc('PT001','DK');


%% 
disp('check_iELVisInstall.m completed. Basic functionality appears to be working.');
