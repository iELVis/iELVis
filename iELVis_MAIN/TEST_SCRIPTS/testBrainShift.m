% This script runs the brain shift correction commands on this page of the
% wiki:
%  http://episurg.pbworks.com/w/page/101851831/Electrode%20Localization

%% Remove existing files
sub='PT001';
eReconDir=fullfile(getFsurfSubDir(),sub,'elec_recon');
eReconFiles{1}='.CT';
eReconFiles{2}='.LEPTO';
eReconFiles{3}='.LEPTOVOX';
eReconFiles{4}='.electrodeNames';
eReconFiles{5}='.INF';
eReconFiles{6}='.PIAL';
eReconFiles{7}='.PIALVOX';
eReconFiles{8}='PostimpLoc.txt';
eReconFiles{9}='.DURAL';
eReconFiles{10}='.DURALVOX';
disp('Running the following commands to clear old files:');
cmnd=sprintf('rm %s',fullfile(eReconDir,'localization_process*'));
disp(cmnd);
[s, w]=unix(cmnd);
for a=1:length(eReconFiles)
    cmnd=sprintf('rm %s',fullfile(eReconDir,[sub eReconFiles{a}]));
    disp(cmnd);
    [s, w]=unix(cmnd);
end


%% Test Yang-Wang brain shift correction code
makeIniLocTxtFile(sub);
yangWangElecPjct(sub);
% Note, 1st grid is 8x8, STP grid is 4 x 5 [1 5 20 16]


%% Remove Files Again
eReconDir=fullfile(getFsurfSubDir(),sub,'elec_recon');
eReconFiles{1}='.CT';
eReconFiles{2}='.LEPTO';
eReconFiles{3}='.LEPTOVOX';
eReconFiles{4}='.electrodeNames';
eReconFiles{5}='.INF';
eReconFiles{6}='.PIAL';
eReconFiles{7}='.PIALVOX';
eReconFiles{8}='PostimpLoc.txt';
eReconFiles{9}='.DURAL';
eReconFiles{10}='.DURALVOX';
disp('Running the following commands to clear old files:');
cmnd=sprintf('rm %s',fullfile(eReconDir,'localization_process*'));
disp(cmnd);
[s, w]=unix(cmnd);
for a=1:length(eReconFiles)
    cmnd=sprintf('rm %s',fullfile(eReconDir,[sub eReconFiles{a}]));
    disp(cmnd);
    [s, w]=unix(cmnd);
end

%% Test Dykstra brain shift correction code
dykstraElecPjct('PT001');


%% Test Plots: Mgrid on Pial Surf
plotMgridOnPial('PT001',1);


%% Test Plot: Mgrid on Slices
cfg=[]; cfg.printFigs=1;
plotMgridOnSlices('PT001',cfg);


%%
disp('Script testBrainShift.m completed successfully.')