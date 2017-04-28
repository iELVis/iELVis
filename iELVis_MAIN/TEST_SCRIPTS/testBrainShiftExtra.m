%% Remove existing files
% Test depths only subject
sub='S13_50_LGM';
eReconDir=fullfile(getFsurfSubDir(),sub,'elec_recon');
eReconFiles{1}='.CT';
eReconFiles{2}='.DURAL';
eReconFiles{3}='.DURALVOX';
eReconFiles{4}='.electrodeNames';
eReconFiles{5}='.INF';
eReconFiles{6}='.PIAL';
eReconFiles{7}='.PIALVOX';
eReconFiles{8}='PostimpLoc.txt';
disp('Running the following commands to clear old files:');
cmnd=sprintf('rm %s',fullfile(eReconDir,'localization_process*'));
disp(cmnd);
[s, w]=unix(cmnd);
for a=1:length(eReconFiles)
    cmnd=sprintf('rm %s',fullfile(eReconDir,[sub eReconFiles{a}]));
    disp(cmnd);
    [s, w]=unix(cmnd);
end


%% Test Dykstra method with depths-only
dykstraElecPjct(sub);

%% Remove files again
eReconDir=fullfile(getFsurfSubDir(),sub,'elec_recon');
eReconFiles{1}='.CT';
eReconFiles{2}='.DURAL';
eReconFiles{3}='.DURALVOX';
eReconFiles{4}='.electrodeNames';
eReconFiles{5}='.INF';
eReconFiles{6}='.PIAL';
eReconFiles{7}='.PIALVOX';
eReconFiles{8}='PostimpLoc.txt';
disp('Running the following commands to clear old files:');
cmnd=sprintf('rm %s',fullfile(eReconDir,'localization_process*'));
disp(cmnd);
[s, w]=unix(cmnd);
for a=1:length(eReconFiles)
    cmnd=sprintf('rm %s',fullfile(eReconDir,[sub eReconFiles{a}]));
    disp(cmnd);
    [s, w]=unix(cmnd);
end


%% Test Yang, Wang method with depths only
makeIniLocTxtFile('S13_50_LGM');
yangWangElecPjct('S13_50_LGM');

%%
disp('Script testBrainShiftExtra.m completed successfully.')
