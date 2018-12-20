%%
iELVisFsurf2BIDS('PT001','~/Desktop/HandMotor',1);

%%
yangWangElecPjct('PT001','~/Desktop/HandMotor',1);

%%
dykstraElecPjct('PT001',0,'~/Desktop/FaceMotor',1);

%% plotPialSurf test
cfg=[];
cfg.view='l';
cfg.figId=1;
cfg.elecCoord='LEPTO';
cfg.title='PT001: iEEG-BIDS';
cfg.bidsDir='/Users/davidgroppe/GIT/iELVis/PRIVATE_FILES/HandMotor'; % Note this differs from path on wiki
cfgOut=plotPialSurf('PT001',cfg);

%% plotMgridOnSlices test
cfg=[]; 
cfg.printFigs='~/Desktop/SLICE_FIGS';
cfg.bidsDir='/Users/davidgroppe/GIT/iELVis/PRIVATE_FILES/HandMotor'; % Note this differs from path on wiki
plotMgridOnSlices('PT001',cfg);