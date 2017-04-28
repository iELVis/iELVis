function infCoor=pial2InfBrain(fsSub,cfg)
%function infCoor=pial2InfBrain(fsSub,cfg)
%
% This function takes the "pial" coordinates (snapped to pial surface) and:
% 1. finds the closest vertex on the pial surface
% 2. maps that vertex to the inflated brain surface
%
% Inputs:
%   fsSub = FreeSurfer subject name
%
% Optional Inputs: passed as fields in a cfg structure
%   elecCoord = N-by-4 numeric array with electrode RAS coordinates. The 
%               fourth column is binary (1=Left hem elec) {default:
%               not used; the function looks into the fsSubect's Freesurfer
%               folder for electrode coordinate files instead}
%   elecNames = cell array of strings with electrode names, corresponding
%               to the rows of elecCoord. {default: not used; the function
%               looks into the fsSubect's Freesurfer folder for electrode
%               name file instead}
%   fsurfsubdir = path to the Freesurfer fsSubect directory. Necessary if
%                 running MATLAB on Windows. {default: taken from shell}
%
% Outputs:
%   infCoords = Electrode coordinates on FreeSurfer inflated pial surface
%
%
% Author: 
% David Groppe
% Mehtalab
% April, 2013
%

% History:
% 2015-6 Made compatible with new Yang, Wang method for brain shift correction: DG

% parse input parameters in cfg structure and set defaults
if  ~isfield(cfg,'elecCoord'),      elecCoord = []; else    elecCoord = cfg.elecCoord;      end
if  ~isfield(cfg,'elecNames'),      elecNames = []; else    elecNames = cfg.elecNames;      end
if  ~isfield(cfg,'fsurfsubdir'),    fs_dir = [];    else    fs_dir = cfg.fsurfsubdir;       end
%if  ~isfield(cfg,'elecHem'),        elecHem = [];   else    elecHem = cfg.elecNames;      end

% Get location of FreeSurfer directories
if isempty(fs_dir)
    fs_dir=getFsurfSubDir();
end
sub_dir=fullfile(fs_dir,fsSub);


%% get electrode coordinates
if isempty(elecCoord) % no electrode coordinates have been passed in the function call:
    % use the original code looking for .PIAL files
    %pialFname=[fs_dir '/' fsSub '/elec_recon/' fsSub '.PIAL'];
    pialFname=fullfile(fs_dir,fsSub,'elec_recon',[fsSub '.PIAL']);
    elecCoordCsv=csv2Cell(pialFname,' ',2);
    nChan=size(elecCoordCsv,1);
    RAS_coor=zeros(nChan,3);
    for csvLoopA=1:nChan,
        for csvLoopB=1:3,
            RAS_coor(csvLoopA,csvLoopB)=str2num(elecCoordCsv{csvLoopA,csvLoopB});
        end
    end
    %elecInfoFname=[fs_dir '/' fsSub '/elec_recon/' fsSub '.electrodeNames'];
    elecInfoFname=fullfile(fs_dir,fsSub,'elec_recon',[fsSub '.electrodeNames']);
    elecInfo=csv2Cell(elecInfoFname,' ',2);
    %labels=elecInfo(:,1);
    leftIds=find(cellfun(@(x) strcmpi(x,'L'),elecInfo(:,3)));
    rightIds=find(cellfun(@(x) strcmpi(x,'R'),elecInfo(:,3)));
else % numeric electrode coordinates have been passed in the function call
    RAS_coor=elecCoord(:,1:3);
    leftIds=find(elecCoord(:,4));
    rightIds=find(~elecCoord(:,4));
    nChan=size(elecCoord,1);
    %     leftIds=find(cellfun(@(x) strcmpi(x,'L'),elecHem));
    %     rightIds=find(cellfun(@(x) strcmpi(x,'R'),elecHem));
    %labels=elecNames;
    elecInfo=[];
end


%% Loop over hemispheres
infCoor=zeros(nChan,3);
for hemLoop=1:2,
    if hemLoop==1,
        hem='l';
        useIds=leftIds;
    else
        hem='r';
        useIds=rightIds;
    end
    if size(useIds,1)>1,
        useIds=useIds';
    end
    
    %% Read Sub Pial Surf
    %fname=[sub_dir '/surf/' hem 'h.pial'];
    fname=fullfile(sub_dir,'surf',[hem 'h.pial']);
    pial=readSurfHelper(fname);
    
    %% Read Sub Inflated Surf
    %fname=[sub_dir '/surf/' hem 'h.inflated'];
    fname=fullfile(sub_dir,'surf',[hem 'h.inflated']);
    inflated=readSurfHelper(fname);
    
    %% Get vertices for electrodes on sub's pial surface
    nPialVert=size(pial.vert,1);
    nUseId=length(useIds);
    if nUseId
        for a=useIds,
            df=pial.vert-repmat(RAS_coor(a,:),nPialVert,1);
            dst=sum(abs(df),2);
            [dummy, sub_vid]=min(dst);
            
            %% Get electrode coordinates on inflated brain
            infCoor(a,:)=inflated.vert(sub_vid,:);
        end
        
        %% Plot to check
        if 0
            figure(10); clf;
            subplot(121);
            plot3(inflated.vert(:,1),inflated.vert(:,2),inflated.vert(:,3),'b.'); hold on;
            plot3(infCoor(useIds,1),infCoor(useIds,2),infCoor(useIds,3),'ro');
            axis square;
            
            subplot(122);
            plot3(pial.vert(:,1),pial.vert(:,2),pial.vert(:,3),'b.'); hold on;
            plot3(RAS_coor(useIds,1),RAS_coor(useIds,2),RAS_coor(useIds,3),'ro');
            axis square;
        end
    end
end

%% Set depth coordinates to NaNs if coordinates read from file
if ~isempty(elecInfo),
    for a=1:size(elecInfo,1),
       if strcmpi(elecInfo{a,2},'D')
          infCoor(a,:)=NaN; 
       end
    end
end
