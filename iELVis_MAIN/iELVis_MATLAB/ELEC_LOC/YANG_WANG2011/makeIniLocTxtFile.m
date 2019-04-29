function makeIniLocTxtFile(fsSub)
%function makeIniLocTxtFile(fsSub)
%
% Input:
%  fsSub - Name of patient's FreeSurfer folder (e.g., fsSub)
%
% Create a text file of elec coordinates (in voxel space) readable by Wang, Yang or Dykstra
% brain shift correctioncode. The file is called *PostimpLoc.txt and is the
% the elec_recon subfolder of the patient's FreeSurfer directory. Note that
% coordinates start at 0 (not 1) and are in the order R->L, S->I, A->P
%
% Author: David M. Groppe
% June, 2015

fsDir=getFsurfSubDir();

% subPath = sprintf('%s/%s',fsDir,fsSub);
% elecReconPath=[subPath '/elec_recon/'];
subPath = fullfile(fsDir,fsSub);
elecReconPath=fullfile(subPath,'elec_recon');
postimpLocFname=fullfile(elecReconPath,[fsSub 'PostimpLoc.txt']);

%% space delimited from mgrid
elecReconDir=fullfile(fsDir,fsSub,'elec_recon');
iLocInfoFname=fullfile(elecReconDir,'iLocElecInfo.tsv');
iLocPairsFname=fullfile(elecReconDir,'iLocElecPairs.tsv');
if exist(iLocInfoFname,'file') && exist(iLocPairsFname,'file')
    % import from iLoc
    [eCoords, elecLabels, elecRgb, elecPairs, elecPresent]=iLoc2Matlab(fsSub);
else
    % import from mgrid
    [eCoords, elecLabels, elecRgb, elecPairs, elecPresent]=mgrid2matlab(fsSub,0);
end
eCoords=eCoords-1; % Make coordinates same as in mgrid file (thus first slice has a coordinate of 0, last has a coordinate of 255)
fprintf('Creating file: %s\n',postimpLocFname);
fid=fopen(postimpLocFname,'w');
for a=1:length(elecLabels)
    erFlag=0;
    if length(elecLabels{a})<4
        erFlag=1;
    elseif ~strcmpi(elecLabels{a}(1),'R') && ~strcmpi(elecLabels{a}(1),'L')
        erFlag=1;
    elseif ~strcmpi(elecLabels{a}(2),'D') && ~strcmpi(elecLabels{a}(2),'S') && ~strcmpi(elecLabels{a}(2),'G')
        erFlag=1;
    end
    if erFlag,
        error('Electrode labels in %s''s mgrid file need to start with RD_, RG_, RS_, LD_, LG_, or LD_',fsSub);
    end
    if strcmpi(elecLabels{a}(1),'R')
        hem='R';
    else
        hem='L';
    end
    if strcmpi(elecLabels{a}(2),'G')
        elecType='G'; %grid
    elseif strcmpi(elecLabels{a}(2),'S')
        elecType='S'; %strip
    else
        elecType='D'; %depth
    end
    
    elecLabels{a}=elecLabels{a}(4:end);
    underIds=find(elecLabels{a}=='_');
    elecNum=str2num(elecLabels{a}((underIds(end)+1):end));
    elecLabels{a}=elecLabels{a}(1:(underIds(end)-1));
    
    fprintf(fid,'%s %d %f %f %f %s %s\n',elecLabels{a},elecNum,eCoords(a,1), ...
        eCoords(a,2),eCoords(a,3),hem,elecType);
end
fclose(fid);
