function plotMgridOnSlices(fsSub,cfg)
% function plotMgridOnSlices(fsSub,cfg)
%
% Creates a figure illustrating the location of each electrode in an mgrid
% file in a sagittal, coronal, and axial slice and indicates which part of
% the brain it is in.
%
% Required Inputs:
%  fsSub - Patient's freesurfer directory name
%
% Optional cfg parameters:
%  mgridFname - mgrid filename and path. If empty, name is assumed to be fsSub.mgrid. 
%  fullTitle  - If 1, the mgrid and mri voxel coordinates are displayed in
%               the figure title along with the electrode name and anatomical 
%               location. {default: 0}
%  markerSize - The size of the dot in each slice used to represent an
%               electrode's location. {default: 30}
%  cntrst     - 0< number <=1 The lower this number, the lower the brightness
%               of the image (i.e., the lower the voxel value corresponding to 
%               white). {default: 0.5}
%  anatOverlay -If 1, color is overlayed on the brain to show FreeSurfer's
%              automatic segmentation of brain areas (neocortex uses 
%              Desikan-Killiany parcellation). Alternatively define the fullpath
%              to another parcellation file. {default: 0}
%  colorLUT    - fullpath to color lookup table if you would like to use 
%                non default colors for your parcellation overlay.
%                {default: FreeSurferColorLUTnoFormat.txt}
%  pauseOn   - If 1, Matlab pa'uses after each figure is made and waits for
%              a keypress. {default: 0}
%  printFigs - If 1, each figure is output to a jpg file in the patient's
%              elec_recon/PICS folder and the figure is closed after the 
%              jpg is created. This is particularly useful for implants with 
%              a large number of depth contacts. {default: 0}
%
%
% Examples:
%  %Specify mgrid file and do NOT print
%  cfg=[];
%  cfg.mgridFname='/Applications/freesurfer/subjects/TWH001/elec_recon/TWH001.mgrid';
%  plotMgridOnSlices('PT001',cfg);
%
%  %Use FreeSurfer file structure and print
%  cfg=[];
%  cfg.printFigs=1;
%  plotMgridOnSlices('PT001',cfg);
%
%
% Author: David M. Groppe
% Feb. 2015
% Feinstein Institute for Medical Research/Univ. of Toronto

% Future work:
% Add option for fsurf anatomy colors?

if ~isfield(cfg,'mgridFname'),    mgridFname=[];    else mgridFname=cfg.mgridFname; end
if ~isfield(cfg,'fullTitle'),     fullTitle=0;      else fullTitle=cfg.fullTitle; end
if ~isfield(cfg,'markerSize'),    markerSize=30;    else markerSize=cfg.markerSize; end
if ~isfield(cfg,'cntrst'),    cntrst=.5;          else cntrst=cfg.cntrst; end
if ~isfield(cfg,'anatOverlay'),    anatOverlay=.5;          else anatOverlay=1; end
if ~isfield(cfg,'colorLUT'),    colorLUT=0;          else colorLUT=cfg.colorLUT; end
if ~isfield(cfg,'pauseOn'),    pauseOn=0;          else pauseOn=cfg.pauseOn; end
if ~isfield(cfg,'printFigs'),    printFigs=0;          else printFigs=cfg.printFigs; end
checkCfg(cfg,'plotMgridOnSlices.m');


% FreeSurfer Subject Directory
fsdir=getFsurfSubDir();

% Load MRI
mriFname=fullfile(fsdir,fsSub,'mri','brainmask.mgz');
if ~exist(mriFname,'file')
   error('File %s not found.',mriFname); 
end
mri=MRIread(mriFname);
%mri.vol is ILA (i.e., S->I, R->L, P->A)
mx=max(max(max(mri.vol)))*cntrst;
mn=min(min(min(mri.vol)));
sVol=size(mri.vol);


% Load mgrid
% if strcmpi(mgridFname,'l') || strcmpi(mgridFname,'r')
%     [elecMatrix, elecLabels, elecRgb]=mgrid2matlab(fsSub,mgridFname);
% else
%     [elecMatrix, elecLabels, elecRgb]=mgrid2matlab(mgridFname); % mgrid coords are LIP
% end
if isempty(mgridFname)
    [elecMatrix, elecLabels, elecRgb]=mgrid2matlab(fsSub);
end
nElec=length(elecLabels);
elecMatrix=round(elecMatrix);
xyz=zeros(size(elecMatrix));
xyz(:,1)=elecMatrix(:,2);
xyz(:,2)=elecMatrix(:,1);
xyz(:,3)=sVol(3)-elecMatrix(:,3);

depthElecs=zeros(nElec,1);
for a=1:nElec,
    if strcmpi(elecLabels{a}(2),'D')
        depthElecs(a)=1;
    end
end


if universalYes(anatOverlay)
    
    % Load segmentation
    if ischar(cfg.anatOverlay)
        segFname = cfg.anatOverlay;
    else
        segFname=fullfile(fsdir,fsSub,'mri','aparc+aseg.mgz');
        if ~exist(mriFname,'file')
            error('File %s not found.',mriFname);
        end
    end
    seg=MRIread(segFname);
    
    % Load segmentation color table
    if universalNo(colorLUT)
        pathstr = fileparts(which('mgrid2matlab'));
        inFile=fullfile(pathstr,'FreeSurferColorLUTnoFormat.txt');
        if ~exist(inFile,'file')
            error('Could not find file %s',inFile);
        end
    elseif exist(colorLUT,'file')
        inFile = colorLUT;
    else
        error('The defined color lookup table was not found');
    end
    fid=fopen(inFile,'r');
    %fid=fopen('/Applications/freesurfer/FreeSurferColorLUTnoFormat.txt','r');
    tbl=textscan(fid,'%d%s%d%d%d%d');
    fclose(fid);
end


for elecId=1:nElec,
    if depthElecs(elecId)
        figId=figure();
        set(figId,'position',[78 551 960 346],'paperpositionmode','auto');
        
        hm=zeros(1,3);
        figure(figId); clf;
        colormap gray;
        %subplot(131);
        wdth=.35;
        wDelt=.33;
        xStart=-.005;
        yStart=.03;
        ht=.9;
        axes('position',[xStart yStart wdth ht]);
        imagesc(squeeze(mri.vol(:,xyz(elecId,2),:)),[mn mx]);
        axis square;
        set(gca,'xdir','reverse');
        hold on;
        
        if universalYes(anatOverlay)
            % Plot segmentation
            for a=1:sVol(1),
                for b=1:sVol(3),
                    if seg.vol(a,xyz(elecId,2),b)
                        segId=find(tbl{1}==seg.vol(a,xyz(elecId,2),b));
                        tempRgb=double([tbl{3}(segId) tbl{4}(segId) tbl{5}(segId)])/255;
                        hM=patch([-.5 .5 .5 -.5]+b,[-.5 -.5 .5 .5]+a,tempRgb);
                        set(hM,'LineStyle','none','FaceAlpha',0.3);
                    end
                end
            end
        end
        
        % Plot electrode
        hm(1)=plot(xyz(elecId,3),xyz(elecId,1),'r.');
        set(hm(1),'color',elecRgb(elecId,:),'markersize',markerSize);
        %find image limits
        mxX=max(squeeze(mri.vol(:,xyz(elecId,2),:)),[],2);
        mxY=max(squeeze(mri.vol(:,xyz(elecId,2),:)),[],1);
        limXa=max(intersect(1:(sVol(3)/2),find(mxX==0)));
        limXb=min(intersect((sVol(3)/2:sVol(3)),find(mxX==0)));
        limYa=max(intersect(1:(sVol(1)/2),find(mxY==0)));
        limYb=min(intersect((sVol(1)/2:sVol(1)),find(mxY==0)));
        %keep image square
        tempMin=min([limXa limYa]);
        tempMax=max([limXb limYb]);
        if tempMin<tempMax, 
            axis([tempMin tempMax tempMin tempMax]); 
        end
        set(gca,'xtick',[],'ytick',[]);
        
        %subplot(132);
        axes('position',[xStart+wDelt yStart wdth ht]);
        imagesc(squeeze(mri.vol(xyz(elecId,1),:,:)),[mn mx]);
        axis square;
        hold on;
        
        if universalYes(anatOverlay)
            % Plot segmentation
            for a=1:sVol(2),
                for b=1:sVol(3),
                    if seg.vol(xyz(elecId,1),a,b)
                        segId=find(tbl{1}==seg.vol(xyz(elecId,1),a,b));
                        tempRgb=double([tbl{3}(segId) tbl{4}(segId) tbl{5}(segId)])/255;
                        hM=patch([-.5 .5 .5 -.5]+b,[-.5 -.5 .5 .5]+a,tempRgb);
                        set(hM,'LineStyle','none','FaceAlpha',0.3);
                    end
                end
            end
        end
        
        hm(2)=plot(xyz(elecId,3),xyz(elecId,2),'r.');
        set(hm(2),'color',elecRgb(elecId,:),'markersize',markerSize);
        %find image limits
        mxX=max(squeeze(mri.vol(xyz(elecId,1),:,:)),[],2);
        mxY=max(squeeze(mri.vol(xyz(elecId,1),:,:)),[],1);
        limXa=max(intersect(1:(sVol(3)/2),find(mxX==0)));
        limXb=min(intersect((sVol(3)/2:sVol(3)),find(mxX==0)));
        limYa=max(intersect(1:(sVol(2)/2),find(mxY==0)));
        limYb=min(intersect((sVol(2)/2:sVol(2)),find(mxY==0)));
        %keep image square
        tempMin=min([limXa limYa]);
        tempMax=max([limXb limYb]);
        if tempMin<tempMax,
            axis([tempMin tempMax tempMin tempMax]);
        end
        set(gca,'xtick',[],'ytick',[],'xdir','reverse');
        
        
        %subplot(133);
        axes('position',[xStart+wDelt*2 yStart wdth ht]);
        imagesc(squeeze(mri.vol(:,:,xyz(elecId,3))),[mn mx]);
        axis square;
        hold on;
        
        if universalYes(anatOverlay)
            % Plot segmentation
            for a=1:sVol(1),
                for b=1:sVol(2),
                    if seg.vol(a,b,xyz(elecId,3))
                        segId=find(tbl{1}==seg.vol(a,b,xyz(elecId,3)));
                        tempRgb=double([tbl{3}(segId) tbl{4}(segId) tbl{5}(segId)])/255;
                        hM=patch([-.5 .5 .5 -.5]+b,[-.5 -.5 .5 .5]+a,tempRgb);
                        set(hM,'LineStyle','none','FaceAlpha',0.3);
                    end
                end
            end
        end
        
        hm(3)=plot(xyz(elecId,2),xyz(elecId,1),'r.');
        set(hm(3),'color',elecRgb(elecId,:),'markersize',markerSize);
        %find image limits
        mxX=max(squeeze(mri.vol(:,:,xyz(elecId,3))),[],2);
        mxY=max(squeeze(mri.vol(:,:,xyz(elecId,3))),[],1);
        limXa=max(intersect(1:(sVol(3)/2),find(mxX==0)));
        limXb=min(intersect((sVol(3)/2:sVol(3)),find(mxX==0)));
        limYa=max(intersect(1:(sVol(2)/2),find(mxY==0)));
        limYb=min(intersect((sVol(2)/2:sVol(2)),find(mxY==0)));
        %keep image square
        tempMin=min([limXa limYa]);
        tempMax=max([limXb limYb]);
        if tempMin<tempMax,
            axis([tempMin tempMax tempMin tempMax]);
        end
        set(gca,'xtick',[],'ytick',[]);
        
        anatLabel=vox2Seg(xyz(elecId,:),fsSub);
  
        % Remove first 3 characters that indicate hemisphere and electrode
        % type
        formattedLabel=elecLabels{elecId}(4:end);
        formattedLabel=rmChar(formattedLabel,'_'); % remove underscore between electrode stem and #
        
        if universalYes(fullTitle)
            ht=textsc2014([formattedLabel '; mgrid coords(' num2str(elecMatrix(elecId,:)-1) '); fsurf coords(' num2str(xyz(elecId,:)) '); ' anatLabel], ...
                'title');
            set(ht,'fontsize',14,'fontweight','bold');
        else
            ht=textsc2014([formattedLabel '; Anatomical Location: ' anatLabel], ...
                'title');
            set(ht,'fontsize',16,'fontweight','bold');
        end
        set(ht,'position',[.5 .97 0]);
        
        if universalYes(printFigs)
            % Make sure PICS directory exists
            erPath=fullfile(fsdir,fsSub,'elec_recon');
            outPath=fullfile(erPath,'PICS');
            if ~exist(outPath,'dir')
                dirSuccess=mkdir(outPath);
                if ~dirSuccess,
                    error('Could not create directory %s',dirSuccess);
                end
            end
            
            drawnow;
            figFname=fullfile(outPath,sprintf('%s_%sSlices',fsSub,elecLabels{elecId}));
            fprintf('Exporting figure to %s\n',figFname);
            %print(figId,figFname,'-depsc');
            print(figId,figFname,'-djpeg');
            pause(1);
            close(figId);
        end
        
        if universalYes(pauseOn)
            fprintf('Paused. Press any key for next electrode.\n');
            pause;
        end
    end
end
fprintf('Done showing all electrodes.\n');