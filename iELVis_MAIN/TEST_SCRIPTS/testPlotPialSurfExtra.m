%% DK atlas with elecs and labels
cfg=[];
%cfg.view='r';
cfg.view='l';
%cfg.view='omni';
cfg.figId=1;
cfg.overlayParcellation='DK';
cfg.showLabels='y';
cfg.title=[];
cfgOut=plotPialSurf('TWH013',cfg);


%% Inflated with elecs and labels
cfg=[];
cfg.view='l';
cfg.figId=2;
cfg.surfType='inflated';
cfg.showLabels='y';
cfg.title=[];
cfgOut=plotPialSurf('PT001',cfg);


%% Plain Omni
cfg=[];
%cfg.view='r';
cfg.view='omni';
cfg.figId=3;
cfg.showLabels='n';
cfg.title=[];
cfgOut=plotPialSurf('TWH013',cfg);
