%% Just negative
cfg=[];
cfg.view='l';
cfg.figId=1;
cfg.pialOverlay='/Applications/freesurfer/subjects/DoBr/fMRI/session_090711/RHandMotor/RHandMotor.feat/reg_surf-lh-DoBr/stats/zstat2.mgh';
cfg.olayColorScale='justneg';
cfg.olayThresh=-1.96;
cfgOut=plotPialSurf('DoBr',cfg);


%% Inflated, just postive
cfg=[];
cfg.view='l';
cfg.figId=1;
cfg.pialOverlay='/Applications/freesurfer/subjects/DoBr/fMRI/session_090711/RHandMotor/RHandMotor.feat/reg_surf-lh-DoBr/stats/zstat2.mgh';
cfg.surfType='inflated';
cfg.olayColorScale='justpos';
cfg.olayThresh=1.96;
cfgOut=plotPialSurf('DoBr',cfg);


%% absmx
cfg=[];
cfg.view='l';
cfg.figId=1;
cfg.surfType='pial';
%cfg.surfType='inflated';
cfg.olayUnits='z';
cfg.pialOverlay='/Applications/freesurfer/subjects/DoBr/fMRI/session_090711/RHandMotor/RHandMotor.feat/reg_surf-lh-DoBr/stats/zstat2.mgh';
cfgOut=plotPialSurf('DoBr',cfg);


%% absmx
cfg=[];
cfg.view='r';
cfg.figId=1;
cfg.surfType='inflated';
cfg.olayUnits='z';
cfg.pialOverlay='/Applications/freesurfer/subjects/DoBr/fMRI/session_090711/RHandMotor/RHandMotor.feat/reg_surf-rh-DoBr/stats/zstat2.mgh';
cfgOut=plotPialSurf('DoBr',cfg);


%% absmx lomni
cfg=[];
cfg.view='lomni';
cfg.figId=2;
cfg.olayUnits='z';
cfg.pialOverlay='/Applications/freesurfer/subjects/DoBr/fMRI/session_090711/RHandMotor/RHandMotor.feat/reg_surf-lh-DoBr/stats/zstat2.mgh';
cfgOut=plotPialSurf('DoBr',cfg);


%% absmx romni
cfg=[];
cfg.view='romni';
cfg.figId=2;
cfg.olayUnits='z';
cfg.pialOverlay='/Applications/freesurfer/subjects/DoBr/fMRI/session_090711/RHandMotor/RHandMotor.feat/reg_surf-rh-DoBr/stats/zstat2.mgh';
cfgOut=plotPialSurf('DoBr',cfg);


%% absmx omni
cfg=[];
cfg.view='omni';
cfg.figId=2;
cfg.olayUnits='z';
cfg.pialOverlay{1}='/Applications/freesurfer/subjects/DoBr/fMRI/session_090711/RHandMotor/RHandMotor.feat/reg_surf-lh-DoBr/stats/zstat2.mgh';
cfg.pialOverlay{2}='/Applications/freesurfer/subjects/DoBr/fMRI/session_090711/RHandMotor/RHandMotor.feat/reg_surf-rh-DoBr/stats/zstat2.mgh';
cfgOut=plotPialSurf('DoBr',cfg);


%% absmx lomni WITH electrode values
elecNames=cell(8,1);
for a=1:8,
    elecNames{a}=sprintf('Gd%d',a+16);
end
elecColors=rand(8,1);

cfg=[];
cfg.view='lomni';
cfg.figId=3;
cfg.olayUnits='z';
cfg.elecColors=elecColors;
cfg.elecUnits='r';
cfg.elecNames=elecNames;
cfg.pialOverlay='/Applications/freesurfer/subjects/DoBr/fMRI/session_090711/RHandMotor/RHandMotor.feat/reg_surf-lh-DoBr/stats/zstat2.mgh';
cfgOut=plotPialSurf('DoBr',cfg);


%% absmx omni WITH electrode values
elecNames=cell(8,1);
for a=1:8,
    elecNames{a}=sprintf('Gd%d',a+16);
end
elecColors=rand(8,1);

cfg=[];
cfg.view='omni';
cfg.figId=3;
cfg.olayUnits='z';
cfg.elecColors=elecColors;
cfg.elecUnits='r';
cfg.elecNames=elecNames;
cfg.pialOverlay{1}='/Applications/freesurfer/subjects/DoBr/fMRI/session_090711/RHandMotor/RHandMotor.feat/reg_surf-lh-DoBr/stats/zstat2.mgh';
cfg.pialOverlay{2}='/Applications/freesurfer/subjects/DoBr/fMRI/session_090711/RHandMotor/RHandMotor.feat/reg_surf-rh-DoBr/stats/zstat2.mgh';
cfgOut=plotPialSurf('DoBr',cfg);


%% absmx l WITH electrode values
elecNames=cell(8,1);
for a=1:8,
    elecNames{a}=sprintf('Gd%d',a+16);
end
elecColors=rand(8,1);

cfg=[];
cfg.view='l';
cfg.figId=3;
cfg.olayUnits='z';
cfg.elecColors=elecColors;
cfg.elecUnits='r';
cfg.elecNames=elecNames;
cfg.pialOverlay='/Applications/freesurfer/subjects/DoBr/fMRI/session_090711/RHandMotor/RHandMotor.feat/reg_surf-lh-DoBr/stats/zstat2.mgh';
cfgOut=plotPialSurf('DoBr',cfg);



%%
disp('testImagingOverlaysExtra.m completed successfully.');