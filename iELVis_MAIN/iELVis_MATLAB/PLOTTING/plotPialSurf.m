% plotPialSurf() - Function for plotting Freesurfer surfaces
%                  with or without colored overlays or electrodes.
%
% Usage:
%  >> cfgOut=plotPialSurf(fsSub,cfg);
%
% Required Input:
%   fsSub - Name of the subject's freesurfer directory (full path not
%          needed)
%
% Optional Inputs:
%   cfg variable with the following possible fields:
%
%    Electrode Options:
%     elecCoord            -If 'n', no electrodes will be rendered in the
%                           figure.  Alternative, you can specify 'LEPTO',
%                           'CT','PIAL', or 'INF' to use the coordinates
%                           with those extensions in the patients elec_recon
%                           folder. *.LEPTO file in patient's
%                           FreeSurfer folder.  Alternatively, you
%                           can pass a 2D matrix of coordinates
%                           instead. The first 3 columns of such a matrix
%                           should be RAS coordinates and the fourth column
%                           is a binary matrix that is 1 if the electrode
%                           is on/in the left hemisphere. {default: 'LEPTO'
%                           for non-inflated brain, 'INF' for inflated brain}
%     elecSize             -Size of electrode markers (disks or spheres).
%                           This also determines thickness of lines connecting
%                           electrodes (if any) and electrode labels (if
%                           shown). {default=8};
%     elecShape            -'marker' or 'sphere': The shape used to
%                           represent electrodes. {default: 'marker'}
%     elecColors           -2D matrix of colors to fill electrodes
%                           (rows=electrodes, columns=RGB values), a column 
%                           vector of values that will be automatically converted
%                           into a color scale, or 'r' to make all red. If
%                           a matrix, the number of rows needs to equal the
%                           number of elecNames (see below). {default: 
%                           all electrodes filled with black}.
%     edgeBlack            -If 'y', electrodes will all have a black
%                           border. Otherwise, border will be same color as
%                           marker. This argument has no effect if
%                           electrodes are represented as spheres. {default:
%                           'y'}
%     elecNames            -Cell array of the names of the electrodes to
%                           show. If elecCoord is a matrix of coordinates,
%                           the number of names needs to equal the number of
%                           rows in the matrix. Otherwise, these names need
%                           to match the names of the electrodes in the
%                           *.electrodeNames file in the subjects Freesurfer
%                           elec_recon folder. If elecCoord specifies a type of
%                           electrode (e.g. LEPTO), then not specifying 
%                           elecNames will cause all electrodes to be plot 
%                           by default.
%     clickElec            -If 'y', clicking on electrodes will reveal
%                           their names in a textbox. Clicking on the box
%                           should make it disapper. Disabled if 'n'. {default: 'y'}
%     ignoreDepthElec      -'y' or 'n': If 'y', depth electrodes will not
%                           be shown. {default: 'y'}
%     pullOut              -Factor via which to project electrodes out from
%                           the center of view. Helpful for when electrodes
%                           sink into the cortical surface. {default: 1}
%     showLabels           -'y' on 'n': If 'y', the name of the first and
%                           last electrode of each strip and each corner of
%                           the 64 chan grid will be shown next to electrode.
%                           {default: 'y'}
%     pairs                -A nx4, nx5, or nx6 cell array specifying n pairs of
%                           electrodes to be connected with lines.
%                           The first two columns indicate which electrodes
%                           are in the pair.
%                           The third column is a 3 element vector indicating
%                           the RGB color of the line that will be drawn
%                           to join the pair of electrodes.
%                           The fourth column is 'l' or 'r' to indicate
%                           which hemisphere the electrodes are on.
%                           The fifth, optional, column is the text that
%                           will appear when the line joining the electrodes
%                           is clicked on.
%                           The sixth, optional, column is the connection
%                           strength for each pair. The maximum value will
%                           correspond to lineWidth; others will be a ratio
%                           of that value.
%                           {default: not used}
%     lineWidth            -Thickness of line connecting pairs of
%                           electrodes. {default: elecSize/3}
%     elecCbar             -'y' or 'n': Plot colorbar next to brain. {default:
%                           'y' if funcfname elecColors argument specified,
%                           'n' otherwise}
%     elecColorScale       -'absmax','minmax', 'justpos', 'justneg',
%                           or numeric vector [minval maxval].  The limits
%                           that define the electrode data color scale.
%                           {default: 'absmax'}.
%     elecUnits             -A string or []. The title of the colorbar
%                           (e.g., 'z-score'). If empty, not title is drawn.
%                           {default: []}
%
%    Surface Options:
%     surfType             -'pial' or 'inflated': Type of Freesurfer surface
%                           to plot. If inflated, gyri and sulci are marked
%                           dark and light grey. {default='pial'};
%     overlayParcellation  -If 'DK', Freesurfer Desikan-Killiany (36 area)
%                           cortical parcellation is plotted on the surface.
%                           If 'D', Freesurfer Destrieux (76 area)
%                           parcellation is used. If 'Y7' or 'Y17', Yeo 
%                           7 or 17-network atlas is used, respectively.
%                           You need to first run createIndivYeoMapping.m
%                           on the individual's data to create Y7 and & Y17 
%                           mappings as they are not produced by default by 
%                           recon-all. You can also input an annotation file
%                           of your choice with fulpath. Be aware that if you 
%                           plot with view = 'omni', you will have to input 
%                           annotation files for the left and right hemisphere 
%                           as a 1x2 cell array, and the filenames have to 
%                           contain 'lh_' and 'rh_'.  {default: not used}
%     parcellationColors   -Optional N-by-3 numeric array with RGB indexes
%                           (0:255) for each of the ROIs in the
%                           overlayParcellation. The colors for each ROI
%                           need to be in the exact same order as that of
%                           the default Freesurfer color table (you can get
%                           that by using the Freesurfer-MATLAB function
%                           [averts,albl,actbl=read_annotation(fname.annot);
%                           the ROIs are listed in actbl.struct_names and
%                           their numeric labels (found in albl) in
%                           actbl.table).
%                           {default: not used; the default Freesurfer
%                           colors for the parcellation are used instead}
%     opaqueness           -[1 to 0] The "alpha" level of the pial surface.
%                           0 means the surface is completely transparent.
%                           1 means that it is completely opaque.
%                           {default: 1}
%     view                 -Angle and lighting with which to view brain.
%                           This also defines which hemisphere to plot.
%                           Options include:
%                             'omni' - 6 views of each hemisphere. When you
%                             use this option with the funcfname option,
%                             funcfname needs to be a cell array of two
%                             filenames. The first specifies the left hem
%                             values and the second the right.
%                             'lomni' - 6 views of left hemisphere
%                             'l' - Left hem lateral
%                             'lm' - Left hem medial
%                             'lo' - Left hem occipital
%                             'lf' - Left hem frontal
%                             'lim' - Left hem inferior-medial
%                             'li' - Left hem inferior
%                             'ls' - Left hem superior
%                             'lsv' - Left hem superior and vertically
%                                     aligned
%                             'liv' - Left hem inferior and vertically
%                                     aligned
%                           Replace 'l' with 'r' to get these views of the
%                           right hemisphere. Alternatively, you can define
%                           everything yourself like so:
%                                   brainView.light=[1 0 0];
%                                   brainView.hem='r';
%                                   brainView.eyes=[45 0]
%                                   cfg.view=brainView
%
%    Neuroimaging Options:
%     pialOverlay          -Filename or cell array of filenames of
%                           functional data to plot.  If plotting both 
%                           hemispheres, use a cell array with the left hem
%                           filename as the 1st element and the right as
%                           the 2nd element. Be sure to include file paths. 
%                           {default: not used}
%     olayColorScale       -'absmax', 'minmax', 'justpos', 'justneg',
%                           or numeric vector (i.e., [minval maxval]).  The 
%                           limits that define the overlay data color scale. 
%                           {default: 'absmax'}
%     olayThresh           -Overlay data with an absolute value less
%                           than this threshold will not be shown (when
%                           showing positive and negative values). If
%                           olayColorScale is 'justpos' then only values
%                           greater than olayThresh will be shown. If
%                           olayColorScale is 'justneg' then only values less
%                           than olayThresh will be shown. {default: not used}
%     olayUnits            -A string or []. The title of the overlay colorbar
%                           (e.g., 'z-score'). If empty, not title is drawn.
%                           {default: []}
%     olayCbar             -'y' or 'n': Plot colorbar next to brain. {default:
%                           'y' if pialOverlay argument specified, 'n' otherwise}
%
%
%    Other Options:
%     axis                 -Handle of axis in which to make plot.
%                           {default: new axis created}
%     figId                -Handle of figure in which to make plot.
%                           {default: new figure created}
%     clearFig             -'y' or 'n'; If 'y', the figured is cleared
%                           before the plot is created. {default: 'y'}
%     backgroundColor      -Standard Matlab color argument (e.g., 'k' or
%                           [.1 .5 .3]). The axis background color.
%                           {default: not used}
%     title                -Title to place on figure. {default: 'y'}
%     fsurfSubDir          -The path to the FreeSurfer subject directory.
%                           Necessary if running MATLAB on Windows.
%                           {default: taken from shell}
%     clearGlobal          -If 'n', some plotting information is left in
%                           global memory. Useful for speeding *omni plots.
%                           {default: 'y'}
%     verbLevel            - An integer specifying the amount of information you want
%                           this function to provide about what it is doing during runtime.
%                            Options are:
%                             0 - quiet, only show errors, warnings, and external function reports
%                             1 - stuff anyone should probably know
%                             2 - stuff you should know the first time you start working
%                                 with a data set {default value}
%                             3 - stuff that might help you debug (show all
%                                 reports)
%
% Example:
% % Plot electrodes on brain with Desikan-Killiany cortical parcellation
% cfg=[];
% cfg.view='l';
% cfg.figId=1;
% cfg.overlayParcellation='DK';
% cfg.showLabels='y';
% cfg.title=[];
% cfgOut=plotPialSurf('PT001',cfg);
%
% % Plot depths with a semi-transparent pial surface
% cfg=[];
% cfg.view='li';
% cfg.ignoreDepthElec='n';
% cfg.opaqueness=.5;
% cfg.title=[];
% cfgOut=plotPialSurf('PT001',cfg);
%
% % Plot electrodes as spheres, color coded to reflect correlation value
% elecNames=cell(6,1);
% for a=1:6,
%     elecNames{a}=sprintf('LO%d',a);
% end
% cfg=[];
% cfg.view='l';
% cfg.figId=1;
% cfg.elecShape='sphere';
% cfg.elecColors=rand(6,1);
% cfg.elecColorScale='minmax';
% cfg.showLabels='n';
% cfg.elecUnits='r';
% cfg.elecNames=elecNames;
% cfg.elecSize=2;
% cfg.title='PT001: Stimulus Correlations';
% cfgOut=plotPialSurf('PT001',cfg);
%
% % Plot electrodes as spheres, strip LO is colored red
% elecNames=cell(6,1);
% for a=1:6,
%     elecNames{a}=sprintf('LO%d',a);
% end
% cfg=[];
% cfg.view='l';
% cfg.figId=1;
% cfg.elecShape='sphere';
% cfg.elecColors=repmat([1, 0, 0],6,1);
% cfg.elecCbar='n';
% cfg.showLabels='n';
% cfg.elecUnits='r';
% cfg.elecNames=elecNames;
% cfg.elecSize=2;
% cfg.title='PT001: Stimulus Correlations';
% cfgOut=plotPialSurf('PT001',cfg);
%
% % Plot bars between electrodes, color coded to reflect bipolar reference correlation value
% pairs=[];
% ct=0;
% for a=1:7,
%     ct=ct+1;
%     pairs{ct,1}=sprintf('Grid%d',a);
%     pairs{ct,2}=sprintf('Grid%d',a+1);
%     pairs{ct,3}=rand(1,3); % RGB val
%     pairs{ct,4}='L';
% end
% for a=1:7,
%     ct=ct+1;
%     pairs{ct,1}=sprintf('Grid%d',a+8);
%     pairs{ct,2}=sprintf('Grid%d',a+1+8);
%     pairs{ct,3}=rand(1,3); % RGB val
%     pairs{ct,4}='L';
% end
% cfg=[];
% cfg.view='l';
% cfg.figId=2;
% cfg.pairs=pairs;
% cfg.showLabels='n';
% cfg.elecUnits='r';
% cfg.title='PT001: Stimulus Correlations';
% cfgOut=plotPialSurf('PT001',cfg);
%
% % Overlay fMRI statistical map from FreeSurfer mgh file
% cfg=[];
% cfg.view='lomni';
% cfg.figId=2;
% cfg.olayUnits='z';
% cfg.pialOverlay='handMotorLH.mgh';
% cfgOut=plotPialSurf('PT001',cfg);
%
%
%  Authors:
%   David M. Groppe, Stephan Bickel, Pierre Mégevand, Andrew Dykstra
%   Laboratory for Multimodal Human Brain Mapping
%   Feinstein Institute for Medical Research
%   Manhasset, New York
%


% HISTORY:
% adykstra/ck 08-2010; get_loc_snap_mgh
% sb 05-2011: do not save graph option
% sb 06-2011: added config structure as input
% dg 03-2011: massively re-written
% dg 5/15/2011: ESM bars between electrodes now also affected by "pullOut"
%  option
% dg 1/3/2013: Updated comments and got rid of underscores in most arguments.
%  Also slightly modified SB's modifications to funcfname argument so that
%  you can pass a vector of values instead of a filename.
% dg 5/13/13: Comments to "pairs" option added (and now pairs can have four
%  columns). click_text changed to clickText3D so that text doesn't
%  disappear into brain.
% dg 9/9/13 changed nominal colors from Destrieux palette to
%  distinguishable_colors.m
% dg 1/2/14 omni, lomni, and romni views added; cfg.units option added
% mm 24/02/14: pairs can know have 5 columns to specify lineWidth for each
%  pairs
% pm x/xx/14 debug passing 2D numeric array for elecCoord
% pm 7/18/14 added optional parcellationColors input
% pm 7/24/14 debug passing 2D numeric array for elecCoord AND overlayParcellation
% pm 7/28/14 debug passing 2D numeric array for elecCoord AND inflated pial surface
% pm 8/12/14 debug ?omni view
% dg 4/15 elecColors can now be a vector of values from which a elecColorScale
% is automatically derived
% dg 6/15 now expects electrode coordinates and names in YangWang method format
% dg 9/15 elecSize now properly modulates sphere
% As of 2016, revision history will be kept on GitHub


%% TO DO
% do auto colormap for electrodes (both jet and rb)?
% Size of electrodes should be relative to brain size?
% Make elecColors and colorbar work for bipolar lines too

function cfgOut=plotPialSurf(fsSub,cfg)

%% Parse parameters
if ~isfield(cfg, 'elecSize'),       elecSize = 8;          else  elecSize = cfg.elecSize;      end
if ~isfield(cfg, 'snap2surf'),      snap2surf = 0;         else  snap2surf = cfg.snap2surf;      end
if ~isfield(cfg, 'surfType'),       surfType = 'pial';     else  surfType = cfg.surfType;     end
if ~isfield(cfg, 'elecCoord'),      elecCoord= 'LEPTO';      else  elecCoord = cfg.elecCoord;       end
if ~isfield(cfg, 'elecColors'),     elecColors= [];        else  elecColors = cfg.elecColors;        end
if ~isfield(cfg, 'elecColorScale'), elecColorScale='absmax';   else elecColorScale=cfg.elecColorScale; end
if ~isfield(cfg, 'elecCbar'),       elecCbar=[];          else elecCbar=cfg.elecCbar; end
if ~isfield(cfg, 'elecUnits'),      elecUnits=[];            else elecUnits=cfg.elecUnits; end
if ~isfield(cfg, 'pullOut'),        pullOut= 1;            else  pullOut = cfg.pullOut; end
if ~isfield(cfg, 'view'),           brainView= 'l';       else  brainView = cfg.view; end
if ~isfield(cfg, 'axis'),           hAx=[];               else  hAx=cfg.axis; end
if ~isfield(cfg, 'overlayParcellation'), overlayParcellation=0;  else  overlayParcellation=cfg.overlayParcellation; end
if ~isfield(cfg, 'parcellationColors'),  parcellationColors= []; else  parcellationColors= cfg.parcellationColors;  end
if ~isfield(cfg, 'figId'),         hFig=[];              else  hFig=cfg.figId; end
if ~isfield(cfg, 'clearFig'),       clearFig=1;            else  clearFig=cfg.clearFig; end
if ~isfield(cfg, 'title'),          surfTitle='default';  else surfTitle=cfg.title; end
if ~isfield(cfg, 'elecNames'),      elecNames=[];    else elecNames=cfg.elecNames; end
if ~isfield(cfg, 'clickElec'),      clickElec='y'; else clickElec=cfg.clickElec; end
if ~isfield(cfg, 'fsurfSubDir'),  fsDir=[];             else fsDir=cfg.fsurfSubDir; end
if ~isfield(cfg, 'verbLevel'),      verbLevel=2;           else verbLevel=cfg.verbLevel; end
if ~isfield(cfg, 'backgroundColor'), backgroundColor=[]; else backgroundColor=cfg.backgroundColor; end
if ~isfield(cfg, 'pairs'), electrode_pairs=[];  else electrode_pairs=cfg.pairs; end
if ~isfield(cfg, 'lineWidth'), lineWidth=[];  else lineWidth=cfg.lineWidth; end
if ~isfield(cfg, 'ignoreDepthElec'), ignoreDepthElec='y'; else ignoreDepthElec=cfg.ignoreDepthElec; end
if ~isfield(cfg, 'edgeBlack'),      edgeBlack='y';         else edgeBlack=cfg.edgeBlack; end
if ~isfield(cfg, 'elecShape'), elecShape='marker';  else elecShape=cfg.elecShape; end
if ~isfield(cfg, 'showLabels'), showLabels=0;  else showLabels=universalYes(cfg.showLabels); end
if ~isfield(cfg, 'opaqueness'),      opaqueness=1;          else opaqueness=cfg.opaqueness; end
if ~isfield(cfg, 'clearGlobal'),    clearGlobal=1;          else clearGlobal=cfg.clearGlobal; end
if ~isfield(cfg, 'pialOverlay'),    pialOverlay=[];         else pialOverlay=cfg.pialOverlay; end 
if ~isfield(cfg, 'olayColorScale'), olayColorScale='absmax';  else olayColorScale=cfg.olayColorScale; end
if ~isfield(cfg, 'olayThresh'), olayThresh=0;  else olayThresh=cfg.olayThresh; end 
if ~isfield(cfg, 'olayCbar'),       olayCbar=[];          else olayCbar=cfg.olayCbar; end 
if ~isfield(cfg, 'olayUnits'),      olayUnits=[];         else olayUnits=cfg.olayUnits; end

global overlayData elecCbarMin elecCbarMax olayCbarMin olayCbarMax; % Needed for ?omni plots

try 
checkCfg(cfg,'plotPialSurf.m');
elecCmapName=[]; % needed for cfgOut
olayCmapName=[]; % needed for cfgOut

% If matrix of elecColors specified make sure an equal number of elecNames
% specfied
if ~isempty(elecColors) && isnumeric(elecColors),
   if size(elecColors,1)~=length(elecNames),
      error('# of elecColors rows, %d, needs to equal # of elecNames, %d.',size(elecColors,1),length(elecNames)); 
   end
end

% If matrix of elecCoord specified make sure an equal number of elecNames
% specfied
if ~strcmpi(elecCoord,'n') && isnumeric(elecCoord),
   if size(elecCoord,1)~=length(elecNames),
      error('# of elecCoord rows, %d, needs to equal # of elecNames, %d.',size(elecCoord,1),length(elecNames)); 
   end
end

if strcmpi(brainView,'omni')
    cfgOut=plotPialOmni(fsSub,cfg);
    return;
elseif strcmpi(brainView,'lomni') || strcmpi(brainView,'romni')
    cfgOut=plotPialHemi(fsSub,cfg);
    return;
end

% FreeSurfer Subject Directory
if isempty(fsDir)
    fsDir=getFsurfSubDir();
end

% Folder with surface files
subFolder=fullfile(fsDir,fsSub);
if ~exist(subFolder,'dir')
   error('FreeSurfer folder for %s not found.',subFolder); 
end
surfacefolder=fullfile(fsDir,fsSub,'surf');

% Get side of brain to show
if ischar(brainView),
    if findstr(brainView,'l')
        side='l';
    else
        side='r';
    end
else
    side=lower(brainView.hem);
    if ~strcmpi(side,'l') && ~strcmpi(side,'r')
        error('cfg.brainView.hem needs to be ''l'' or ''r''.');
    end
end

% Set/check electrode colorbar plotting options
if isempty(elecCbar)
    if ~isempty(elecColors) && ~ischar(elecColors)
        elecCbar='y';
    else
        elecCbar='n';
    end
end
if isempty(olayCbar)
    if ~isempty(pialOverlay)
        olayCbar='y';
    else
        olayCbar='n';
    end
end

%% %%%%%%%%%% Start Main Function %%%%%%%
verbReport('**** PLOTTING CORTICAL SURFACE WITH "plotPialSurf.m" ****', ...
    2,verbLevel);


%% MAKE FIGURE/AXIS
if ~isempty(hFig),
    figure(hFig);
else
    hFig=figure;
end
if universalYes(clearFig),
    clf;
end
if ~isempty(backgroundColor)
    set(hFig,'color',backgroundColor);
end
if ~isempty(hAx),
    axes(hAx);
else
    hAx=gca;
end


%% If plotting on inflated surface, load curvature values so that sulci and
% gyri can be seen
if strcmpi(surfType,'inflated')
    if side == 'r'
        curv = read_curv(fullfile(surfacefolder,'rh.curv'));
    else
        curv = read_curv(fullfile(surfacefolder,'lh.curv'));
    end
    curvMap=zeros(length(curv),3);
    pcurvIds=find(curv>=0);
    curvMap(pcurvIds,:)=repmat([1 1 1]*.3,length(pcurvIds),1);
    ncurvIds=find(curv<0);
    curvMap(ncurvIds,:)=repmat([1 1 1]*.7,length(ncurvIds),1);
end

%% Initialize pial surface coloration
% Color gyri and sulci different shades of grey
if side == 'r'
    curv = read_curv(fullfile(surfacefolder,'rh.curv'));
else
    curv = read_curv(fullfile(surfacefolder,'lh.curv'));
end
if strcmpi(surfType,'inflated')
    overlayDataTemp=zeros(length(curv),3);
    pcurvIds=find(curv>=0);
    overlayDataTemp(pcurvIds,:)=repmat([1 1 1]*.3,length(pcurvIds),1);
    ncurvIds=find(curv<0);
    overlayDataTemp(ncurvIds,:)=repmat([1 1 1]*.7,length(ncurvIds),1);
else
    overlayDataTemp=ones(length(curv),3)*.7;
end


if ~isempty(overlayData) || ~isempty(pialOverlay)
    % Pial Surface Overlay (e.g., fMRI statistical maps)
    if  ~isempty(pialOverlay)
        [pathstr, name, ext]=fileparts(pialOverlay);
        if strcmpi(ext,'.mgh')
            % FreeSurfer formatted file
            mgh = MRIread(pialOverlay);
            overlayData=mgh.vol;
        else
            % mat file that needs to contain an mgh variable
            if ~exist(pialOverlay,'file')
                error('File %s not found.',pialOverlay);
            end
            load(pialOverlay,'overlayData');
            if ~exist('overlayData','var')
                error('File %s does NOT contain a variable called overlayData.',pialOverlay);
            end
        end
    end
    if isvector(overlayData)
        %overlayData is a vector of values that needs to be converted to
        %RGB
        olayDataVec=overlayData;
        [overlayData, oLayLimits, olayCmapName]=vals2Colormap(olayDataVec,olayColorScale);
        olayCbarMin=oLayLimits(1);
        olayCbarMax=oLayLimits(2);
        if strcmpi(olayColorScale,'justpos')
            maskIds=find(olayDataVec<=olayThresh);
            overlayData(maskIds,:)=overlayDataTemp(maskIds,:); % make subthreshold values grey
            %overlayData(maskIds,:)=repmat([.7 .7 .7],length(maskIds),1); % make subthreshold values grey
        elseif strcmpi(olayColorScale,'justneg')
            maskIds=find(olayDataVec>=olayThresh);
            overlayData(maskIds,:)=overlayDataTemp(maskIds,:); % make superthreshold values grey
            %overlayData(maskIds,:)=repmat([.7 .7 .7],length(maskIds),1); % make superthreshold values grey
        elseif olayThresh~=0
            maskIds=find(abs(olayDataVec)<=olayThresh);
            overlayData(maskIds,:)=overlayDataTemp(maskIds,:); % make abs subthreshold values grey
            %overlayData(maskIds,:)=repmat([.7 .7 .7],length(maskIds),1); % make abs subthreshold values grey
        end
        clear olayDataVec
    end
else
    overlayData=overlayDataTemp;
end
clear overlayDataTemp


%% READ SURFACE
global cort %speeds up omni a tiny bit
if isempty(cort)
    if side == 'r'
        [cort.vert cort.tri]=read_surf(fullfile(surfacefolder,['rh.' surfType]));
    else
        [cort.vert cort.tri]=read_surf(fullfile(surfacefolder,['lh.' surfType]));
    end
    if min(min(cort.tri))<1
        cort.tri=cort.tri+1; %sometimes this is needed sometimes not. no comprendo. DG
    end
end


%% PLOT SURFACE
tripatchDG(cort,hFig,overlayData); %this plots the brain
rotate3d off;


%% If specified, overlay cortical parcellation
if overlayParcellation,
    labelFolder=fullfile(fsDir,fsSub,'label');
    if ~isempty(labelFolder) && (labelFolder(end)~='/')
        labelFolder=[labelFolder '/'];
    end
    if strcmpi(overlayParcellation,'DK')
        annotFname=fullfile(labelFolder,[side 'h.aparc.annot']); %Desikan-Killiany 36 area atlas
        [averts,albl,actbl]=read_annotation(annotFname);
        actbl.table(1,1:3)=255*[1 1 1]*.7; %make medial wall the same shade of grey as functional plots
    elseif strcmpi(overlayParcellation,'D')
        annotFname=fullfile(labelFolder,[side 'h.aparc.a2009s.annot']); %Destrieux 76 area atlas
        [averts,albl,actbl]=read_annotation(annotFname);
        actbl.table(43,1:3)=255*[1 1 1]*.7; %make medial wall the same shade of grey as functional plots
    elseif strcmpi(overlayParcellation,'Y7')
        if strcmpi(fsSub,'fsaverage')
            annotFname=fullfile(labelFolder,[side 'h.Yeo2011_7Networks_N1000.annot']); % Yeo et al. 2011
            [averts,albl,actbl]=read_annotation(annotFname);
        else
            annotFname=fullfile(labelFolder,[side 'h_Yeo2011_7Networks_N1000.mat']); % Yeo et al. 2011
            load(annotFname);
            albl=label;
            actbl=colortable;
            clear colortable label vertices
            actbl.table(1,1:3)=255*[1 1 1]*.7; %make medial wall the same shade of grey as functional plots
        end
    elseif strcmpi(overlayParcellation,'Y17')
        if strcmpi(fsSub,'fsaverage')
            annotFname=fullfile(labelFolder,[side 'h.Yeo2011_17Networks_N1000.annot']); % Yeo et al. 2011
            [averts,albl,actbl]=read_annotation(annotFname);
        else
            annotFname=fullfile(labelFolder,[side 'h_Yeo2011_17Networks_N1000.mat']); % Yeo et al. 2011
            load(annotFname);
            albl=label;
            actbl=colortable;
            clear colortable label vertices
            actbl.table(1,1:3)=255*[1 1 1]*.7; %make medial wall the same shade of grey as functional plots
        end
    elseif exist(overlayParcellation,'file')
        [averts,albl,actbl]=read_annotation(overlayParcellation);
        %  actbl.table(43,1:3)=255*[1 1 1]*.7; %make medial wall the same shade of grey as functional plots
    else
        error('overlayParcellation argument needs to take a value of ''D'',''DK'',''Y7'', ''Y17'', or fullpath to annotation file.');
    end
    if ~isempty(parcellationColors)
        if size(parcellationColors,1)~=size(actbl.table,1)
            error('plotPialSurf:colors_parcellation_size1','parcellationColors argument needs to have \n the same number of rows as the number of ROIs \n in the parcellation. For %s, %d.',overlayParcellation,size(actbl.table,1));
        end
        if size(parcellationColors,2)~=3
            error('plotPialSurf:colors_parcellation_size2','parcellationColors must be an N-by-3 array.');
        end
        actbl.table(:,1:3)=parcellationColors;
    end
    clear averts;
    [~,locTable]=ismember(albl,actbl.table(:,5));
    locTable(locTable==0)=1; % for the case where the label for the vertex says 0
    fvcdat=actbl.table(locTable,1:3)./255; %scale RGB values to 0-1
    clear locTable;
    hTsurf=trisurf(cort.tri, cort.vert(:, 1), cort.vert(:, 2), cort.vert(:, 3),...
        'FaceVertexCData', fvcdat,'FaceColor', 'interp','FaceAlpha',1);
    if ~universalNo(elecCoord),
        cfg_elec2parc=[];
        cfg_elec2parc.fsurfSubDir=fsDir;
        if isnumeric(elecCoord)
            cfg_elec2parc.elecCoord=elecCoord;
            cfg_elec2parc.elecNames=cfg.elecNames;
        end
        elecAssign=[];
    end
else
    elecAssign=[];
end


%% Set Lighting & View
shading interp; lighting gouraud; material dull; axis off, hold on
if ischar(brainView)
    switch brainView
        case 'r'
            l=light('Position',[1 0 0]);
            view(90,0)
        case 'rm'
            l=light('Position',[-1 0 0]);
            view(270,0)
        case 'rim'
            l=light('Position',[-1 0 0]);
            view(270,-45)
        case 'ri'
            l=light('Position',[0 0 -1]);
            view(90,-90)
        case 'ro'
            l=light('Position',[0 -1 0]);
            view(0,0)
        case 'lo'
            l=light('Position',[0 -1 0]);
            view(0,0)
        case 'rf'
            l=light('Position',[0 1 0]);
            view(180,0)
        case 'lf'
            l=light('Position',[0 1 0]);
            view(180,0)
        case 'rs'
            l=light('Position',[0 0 1]);
            view(90,90);
        case 'rsv' %superior & vertically aligned
            l=light('Position',[0 0 1]);
            view(0,90);
        case 'l'
            l=light('Position',[-1 0 0]);
            view(270,0);
        case 'lm'
            l=light('Position',[1 0 0]);
            view(90,0);
        case 'li'
            l=light('Position',[0 0 -1]);
            view(90,-90);
        case 'lim'
            l=light('Position',[-1 0 0]);
            view(270,-45);
        case 'ls'
            l=light('Position',[0 0 1]);
            view(270,90);
        case 'lsv' %superior & vertically aligned
            l=light('Position',[0 0 1]);
            view(0,90);
        case 'liv' %inferior & vertically aligned
            l=light('Position',[0 0 -1]);
            view(0,-90);
        case 'riv' %inferior & vertically aligned
            l=light('Position',[0 0 -1]);
            view(0,-90);
    end
    clear l
else
    light('Position',brainView.light);
    view(brainView.eyes)
end
alpha(opaqueness);


%% PLOT ELECTRODES (optional)
if universalNo(elecCoord)
    verbReport('...not plotting electrodes',2,verbLevel);
    h_elec=[];
else
    [showElecCoords, showElecNames, h_elec, elecCbarMin, elecCbarMax]=plotElecs(elecCoord, ...
        surfType,fsDir,fsSub,side,ignoreDepthElec,pullOut,elecColors,elecColorScale, ...
        elecShape,elecSize,showLabels,clickElec,elecAssign,edgeBlack,elecNames, ...
        elecCbar);
end


%% Add Title
if strcmpi(surfTitle,'default'),
    if ischar(brainView)
        surfTitle=[fsSub '; ' brainView '; '];
    else
        surfTitle=[fsSub '; ' brainView.hem '; '];
    end
    title(surfTitle,'fontsize',20);
elseif ~isempty(surfTitle)
    title(surfTitle,'fontsize',20);
end


%% Colorbar
hElecCbar=[]; % needs to be declared for cfgOut even if colorbar not drawn
hOlayCbar=[]; % needs to be declared for cfgOut even if colorbar not drawn
if universalYes(elecCbar) && universalYes(olayCbar)
    % Plot both electrode AND olay colorbar
    if isempty(elecCbarMin) || isempty(elecCbarMax)
        fprintf('elecCbarMin or elecCbarMax are empty. Cannot draw colorbar.\n');
    else
        pos=[0.9158 0.1100 0.0310 0.8150];
        cbarFontSize=12;
        cbarDGplus(pos,[elecCbarMin elecCbarMax],elecCmapName,5,elecUnits,'top',cbarFontSize);
    end
    
    pos=[0.1 0.05 0.8150 0.0310];
    cbarDGplus(pos,[olayCbarMin olayCbarMax],olayCmapName,5,olayUnits,'top',cbarFontSize);
elseif universalYes(elecCbar)
    % Plot electrode colorbar only
    if isempty(elecCbarMin) || isempty(elecCbarMax)
        fprintf('elecCbarMin or elecCbarMax are empty. Cannot draw colorbar.\n');
    else
        % Plot electrode colorbar only
        pos=[0.9158 0.1100 0.0310 0.8150];
        cbarDGplus(pos,[elecCbarMin elecCbarMax],elecCmapName,5,elecUnits);
        
        if isequal(get(hFig,'color'),[0 0 0]);
            %If background of figure is black, make colorbar text white
            set(hElecCbar,'xcolor','w'); % fix so that box isn't white? ??
            set(hElecCbar,'ycolor','w');
        end
    end
elseif universalYes(olayCbar)
    % Plot pial surface overlay colorbar only
    pos=[0.9158 0.1100 0.0310 0.8150];
    cbarDGplus(pos,[olayCbarMin olayCbarMax],olayCmapName,5,olayUnits);

    if isequal(get(hFig,'color'),[0 0 0]);
        %If background of figure is black, make colorbar text white
        set(hOlayCbar,'xcolor','w'); % fix so that box isn't white? ??
        set(hOlayCbar,'ycolor','w');
    end
end


%% COLLECT CONFIG OUTPUT
cfgOut.subject=fsSub;
cfgOut.view=brainView;
cfgOut.elecSize=elecSize;
cfgOut.elecHandles=h_elec;
cfgOut.surfType=surfType;
cfgOut.hElecCbar=hElecCbar;
cfgOut.hOlayCbar=hOlayCbar;
cfgOut.hBrain=hAx;
cfgOut.elecCmapName=elecCmapName;
cfgOut.olayCmapName=olayCmapName;
if exist('cfg','var'), cfgOut.cfg=cfg; end
if exist('showElecCoords','var'), 
    cfgOut.electrodeCoords=showElecCoords; 
    cfgOut.electrodeNames=showElecNames;
end 
if exist('elecCbarMin','var'), cfgOut.elecCbarLimits=[elecCbarMin elecCbarMax]; end
if exist('olayCbarMin','var'), cfgOut.olayCbarLimits=[olayCbarMin olayCbarMax]; end

if universalYes(clearGlobal)
    clear global elecCbarMin elecCbarMax olayCbarMin olayCbarMax cort overlayData;
end

catch err
    disp(err.identifier);
    disp(err.message);
    for errLoop=1:length(err.stack),
        disp(err.stack(errLoop).file); 
        fprintf('Line: %d\n',err.stack(errLoop).line); 
    end
    % Delete global variables if function crashes to prevent them from
    % being automatically used the next time plotPialSurf is called.
    clear global elecCbarMin elecCbarMax olayCbarMin olayCbarMax cort overlayData;
end

%%%% END OF MAIN FUNCTION %%%%


%% HELPER FUNCTIONS

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% subfunction plotPialOmni
function sub_cfg_out=plotPialOmni(fsSub,cfg)

if ~isfield(cfg, 'figId'),         hFig=[];            else  hFig=cfg.figId; end
if ~isfield(cfg, 'olayThresh'),       olayThresh=[];          else  olayThresh = cfg.olayThresh; end
if ~isfield(cfg, 'figId'),         hFig=[];              else  hFig=cfg.figId; end
if ~isfield(cfg, 'fsurfSubDir'),   fsDir=[];             else fsDir=cfg.fsurfSubDir; end
if ~isfield(cfg, 'elecCoord'),      elecCoord='LEPTO';      else  elecCoord = cfg.elecCoord;       end
if ~isfield(cfg, 'elecSize'),       elecSize=8;          else  elecSize = cfg.elecSize;      end
if ~isfield(cfg, 'elecColors'),     elecColors=[];        else  elecColors = cfg.elecColors;        end
if ~isfield(cfg, 'elecColorScale'),   elecColorScale='absmax';   else elecColorScale=cfg.elecColorScale; end
if ~isfield(cfg, 'olayColorScale'),   olayColorScale='absmax';   else olayColorScale=cfg.olayColorScale; end
if ~isfield(cfg, 'elecUnits'),     elecUnits=[];   else elecUnits=cfg.elecUnits; end
if ~isfield(cfg, 'olayUnits'),      olayUnits=[];         else olayUnits=cfg.olayUnits; end 
if ~isfield(cfg, 'showLabels'),         showLabels='y';            else  showLabels=cfg.showLabels; end
if ~isfield(cfg, 'elecCbar'),     elecCbar=[];   else elecCbar=cfg.elecCbar; end
if ~isfield(cfg, 'olayCbar'),     olayCbar=[];   else elecCbar=cfg.olayCbar; end
if ~isfield(cfg, 'verbLevel'),     verbLevel=0;        else  verbLevel = cfg.verbLevel;        end
if ~isfield(cfg, 'pialOverlay'),    pialOverlay=[];        else pialOverlay=cfg.pialOverlay; end 

if isempty(fsDir)
    fsDir=getFsurfSubDir();
end

clear global elecCbarMin elecCbarMax olayCbarMin olayCbarMax cort;

% Optional electrode color bar variables
if ~isempty(elecColors),
    if isempty(elecCbar)
        elecCbar='y';
    end
end
elecCmapName=[];
if isnumeric(elecColorScale)
    elecUsedLimits=elecColorScale;
else
    elecUsedLimits=[];
end

% Optional pial overlay color bar variables
olayCmapName=[];
olayUsedLimits=[];
if ~isempty(pialOverlay),
    if isempty(olayCbar)
        olayCbar='y';
    end
    if length(pialOverlay)~=2
        error('When you use the "omni" view option and "pialOverlay," pialOverlay needs two filenanes (cfg.pialOverlay{1}=left hem, cfg.pialOverlay{2}=right hem).');
    end
    % Load surf files to get color scale limits
    for a=1:2,
        if ~exist(pialOverlay{a},'file')
            error('File not found: %s',pialOverlay{a});
        end
        dot_id=find(pialOverlay{a}=='.');
        extnsn=pialOverlay{a}(dot_id+1:end);
        if strcmpi(extnsn,'mat')
            load(pialOverlay{a});
        else
            mgh = MRIread(pialOverlay{a});
        end
        switch lower(olayColorScale)
            case 'absmax'
                abs_mx(a)=max(abs(mgh.vol));
            case 'justpos'
                func_mx(a)=max(mgh.vol);
                olayCmapName='autumn';
            case 'justneg'
                func_mn(a)=min(mgh.vol);
                olayCmapName='winter';
            case 'minmax'
                func_mn(a)=min(mgh.vol);
                func_mx(a)=max(mgh.vol);
            otherwise
                %'nominal',
                error('Invalid value of cfg.olayColorScale');
        end
    end
        
    switch lower(olayColorScale)
        case 'absmax'
            olayUsedLimits=[-1 1]*max(abs_mx);
        case 'justpos'
            olayUsedLimits=[olayThresh max(func_mx)];
        case 'justneg'
            olayUsedLimits=[min(func_mn) olayThresh];
        case 'minmax'
            olayUsedLimits=[min(func_mn) max(func_mx)];
    end
end

% Create figure
if isempty(hFig),
    hFig=figure; clf;
else
    figure(hFig); clf;
end
set(hFig,'MenuBar','none','position',[100 190 1000 600],'paperpositionmode','auto');

% Figure out which hemisphere has electrodes
if isnumeric(elecCoord),
    leftCoverage=sum(elecCoord(:,4))>0;
    rightCoverage=sum(~elecCoord(:,4))>0;
else
    % Grab electrode info from subject dir
    elecInfoFname=fullfile(fsDir,fsSub,'elec_recon',[fsSub '.electrodeNames']);
    elecInfo=csv2Cell(elecInfoFname,' ',2);
    leftCoverage=~isempty(findStrInCell('L',elecInfo(:,3)));
    rightCoverage=~isempty(findStrInCell('R',elecInfo(:,3)));
end

% Loop over hemispheres
for h=1:2,
    for v=1:6,
        ax_loc=[0 0 0 0];
        if h==1,
            bview='l';
        else
            bview='r';
        end
        if v==2,
            bview=[bview 'm'];
        elseif v==3,
            bview=[bview 'f'];
        elseif v==4,
            bview=[bview 'o'];
        elseif v==5,
            bview=[bview 'sv'];
        elseif v==6,
            bview=[bview 'iv'];
        end
        switch (h-1)*6+v,
            case 1 %LH lateral
                ax_loc=[0 .67 .4 .3];
            case 2 %LH medial
                ax_loc=[0 .34 .4 .3];
            case 3 %LH frontal
                ax_loc=[.155 .02 .2 .3];
            case 4 %LH occiptal
                ax_loc=[1-.15-.2 .02 .2 .3];
            case 5 %LH superior
                ax_loc=[1-.455-.2 .55 .2 .4];
            case 6 %LH inferior
                ax_loc=[1-.455-.2 .14 .2 .4];
                
            case 7 %RH lateral
                ax_loc=[1-.4 .67 .4 .3];
            case 8 %RH medial
                ax_loc=[1-.4 .34 .4 .3];
            case 9 %RH frontal
                ax_loc=[.045 .02 .2 .3];
            case 10 %RH occipital
                ax_loc=[1-.035-.2 .02 .2 .3];
            case 11 %RH superior
                ax_loc=[.455 .55 .2 .4];
            case 12 %RH inferior
                ax_loc=[.455 .14 .2 .4];
        end
        hAx=axes('position',ax_loc);
        sub_cfg=cfg;
        sub_cfg.view=bview;
        
        sub_cfg.title=[];
        if bview(1)=='l'
            if ~leftCoverage || (ischar(elecCoord) && universalNo(elecCoord)), 
                sub_cfg.elecCoord='n';
            else
                sub_cfg.elecCoord=elecCoord;
                if ~isfield(sub_cfg,'elecSize')
                    sub_cfg.elecSize=6;
                end
                sub_cfg.showLabels=showLabels;
            end
        else
            if ~rightCoverage || (ischar(elecCoord) && universalNo(elecCoord)),
                sub_cfg.elecCoord='n';
            else
                sub_cfg.elecCoord=elecCoord;
                if ~isfield(sub_cfg,'elecSize')
                    sub_cfg.elecSize=6;
                end
                sub_cfg.showLabels=showLabels;
            end
        end
        
        % If plotting annotations from files input by user, choose hemisphere
        if isfield(sub_cfg,'overlayParcellation')
            if iscell(sub_cfg.overlayParcellation) && length(sub_cfg.overlayParcellation)==2
                
                if ~isempty(strfind(sub_cfg.overlayParcellation{1},'lh_'))  && ~isempty(strfind(sub_cfg.overlayParcellation{2},'rh_'))
                    lh_annot = sub_cfg.overlayParcellation{1};
                    rh_annot = sub_cfg.overlayParcellation{2};
                elseif ~isempty(strfind(sub_cfg.overlayParcellation{1},'rh_'))  && ~isempty(strfind(sub_cfg.overlayParcellation{2},'lh_'))
                    lh_annot = sub_cfg.overlayParcellation{2};
                    rh_annot = sub_cfg.overlayParcellation{1};
                else
                    error(['If you define parcellation files to overlay with omni view,'...
                        'cfg.overlayParcellation has to be a cell array with the left and '...
                        'right hemisphere annotation file. The file names (with fullpath) have to include either ''lh_'' or ''rh_''']);
                end
                
                if h==1
                    sub_cfg.overlayParcellation =  lh_annot;
                elseif h==2
                    sub_cfg.overlayParcellation =  rh_annot;
                end
            end
        end
        
        sub_cfg.figId=hFig;
        sub_cfg.axis=hAx;
        sub_cfg.verbLevel=verbLevel;
        if ~isempty(pialOverlay),
            sub_cfg.pialOverlay=pialOverlay{h};
        end
        sub_cfg.clearFig='n';
        sub_cfg.elecCbar='n';
        sub_cfg.olayCbar='n';
        sub_cfg.olayColorScale=olayUsedLimits;
        if v==6
            sub_cfg.clearGlobal=1; %last view for this hem, clear overlay data from global memory
        else
            sub_cfg.clearGlobal=0;
        end
        sub_cfg_out=plotPialSurf(fsSub,sub_cfg);
        
        % Get electrode colormap limits
        
        if isempty(elecUsedLimits)
            if isfield(sub_cfg_out,'elecCbarLimits')
                elecUsedLimits=sub_cfg_out.elecCbarLimits;
            end
        else
            if ~isempty(sub_cfg_out.elecCbarLimits)
                if elecUsedLimits(2)<sub_cfg_out.elecCbarLimits(2)
                    elecUsedLimits(2)=sub_cfg_out.elecCbarLimits(2);
                end
                if elecUsedLimits(1)>sub_cfg_out.elecCbarLimits(1)
                    elecUsedLimits(1)=sub_cfg_out.elecCbarLimits(1);
                end
            end
        end
        
        if isempty(elecCmapName) && isfield(sub_cfg_out,'elecCmapName')
            elecCmapName=sub_cfg_out.elecCmapName;
        end
    end
end


%% DRAW COLORBAR(S)
if universalYes(elecCbar) && universalYes(olayCbar),
    % Colorbar for electrodes
    pos=[.4 .09 .2 .01];
    cbarDGplus(pos,elecUsedLimits,elecCmapName,5,elecUnits,'right');

    % Colorbar for pial surface overlay (e.g., neuroimaging)
    pos=[.4 .04 .2 .01];
    cbarDGplus(pos,olayUsedLimits,olayCmapName,5,olayUnits,'right');
elseif universalYes(elecCbar),
    % Colorbar for electrodes
    pos=[.4 .06 .2 .03];
    cbarDGplus(pos,elecUsedLimits,elecCmapName,5,elecUnits);
elseif universalYes(olayCbar)
    % Colorbar for pial surface overlay (e.g., neuroimaging)
    pos=[.4 .06 .2 .03];
    cbarDGplus(pos,olayUsedLimits,olayCmapName,5,olayUnits);
end

% Title
if isfield(cfg,'title')
    if ~isempty(cfg.title)
        % Overall Fig Title
        ht=textsc2014(cfg.title,'title');
        set(ht,'fontweight','bold','fontsize',20,'position',[0.5 0.975]);
    end
end
drawnow;


%% subfunction plotPialHemi
function sub_cfg_out=plotPialHemi(fsSub,cfg)

if ~isfield(cfg, 'usemask'),       usemask=[];          else usemask=cfg.usemask; end
if ~isfield(cfg, 'figId'),         hFig=[];            else  hFig=cfg.figId; end
if ~isfield(cfg, 'fsurfSubDir'),   fsDir=[];             else fsDir=cfg.fsurfSubDir; end
if ~isfield(cfg, 'elecCoord'),      elecCoord= 'LEPTO';      else  elecCoord = cfg.elecCoord;       end
if ~isfield(cfg, 'elecSize'),       elecSize = 8;          else  elecSize = cfg.elecSize;      end
if ~isfield(cfg, 'elecColors'),     elecColors= [];        else  elecColors = cfg.elecColors;        end
if ~isfield(cfg, 'elecColorScale'),   elecColorScale=[];   else elecColorScale=cfg.elecColorScale; end
if ~isfield(cfg, 'olayColorScale'),   olayColorScale=[];   else olayColorScale=cfg.olayColorScale; end
if ~isfield(cfg, 'elecUnits'),     elecUnits=[];   else elecUnits=cfg.elecUnits; end
if ~isfield(cfg, 'olayUnits'),      olayUnits=[];         else olayUnits=cfg.olayUnits; end 
if ~isfield(cfg, 'elecCbar'),     elecCbar=[];   else elecCbar=cfg.elecCbar; end
if ~isfield(cfg, 'olayCbar'),     olayCbar=[];   else elecCbar=cfg.olayCbar; end
if ~isfield(cfg, 'verbLevel'),     verbLevel=0;        else  verbLevel = cfg.verbLevel;        end
if ~isfield(cfg, 'pialOverlay'),    pialOverlay=[];        else pialOverlay=cfg.pialOverlay; end 

if ~isempty(elecColors),
    if isempty(elecCbar)
        elecCbar='y';
    end
end
elecCmapName=[];
if ~isempty(pialOverlay),
    if isempty(olayCbar)
        olayCbar='y';
    end
end
olayCmapName=[];

clear global elecCbarMin elecCbarMax olayCbarMin olayCbarMax cort;

hem=cfg.view(1);

if isempty(hFig),
    hFig=figure; clf;
else
    figure(hFig); clf;
end
set(hFig,'MenuBar','none','position',[100 190 800 500],'paperpositionmode','auto');

%% Viewpoints
elecUsedLimits=[];
olayUsedLimits=[];
for v=1:6, %will run 1-6
    ax_loc=[0 0 0 0];
    bview=hem;
    if v==2,
        bview=[bview 'm'];
    elseif v==3,
        if bview(1)=='r',
            bview=[bview 'f'];
        else
            bview=[bview 'o'];
        end
    elseif v==4,
        if bview(1)=='r'
            bview=[bview 'o'];
        else
            bview=[bview 'f'];
        end
    elseif v==5,
        bview=[bview 's'];
    elseif v==6,
        bview=[bview 'i'];
    end
    switch v,
        case 1 % lateral
            ax_loc=[-.03 .52 .55 .45];
        case 2 % medial
            ax_loc=[-.03 .05 .55 .45];
        case 3 % occipital
            ax_loc=[.305 .55 .55 .41];
            %ax_loc=[.29 .55 .55 .41];
        case 4 % frontal
            ax_loc=[.56 .55 .44 .41];
            %ax_loc=[.56 .55 .44 .41];
        case 5 % superior
            ax_loc=[.41 .05 .54 .23];
        case 6 % inferior
            ax_loc=[.41 .31 .54 .23];
    end
    hAx=axes('position',ax_loc);
    
    sub_cfg=cfg;
    sub_cfg.view=bview;
    sub_cfg.title=[];
    sub_cfg.figId=hFig;
    sub_cfg.axis=hAx;
    sub_cfg.verbLevel=verbLevel;
    sub_cfg.clearFig='n';
    sub_cfg.elecCbar='n';
    sub_cfg.olayCbar='n';
    if v==6
        sub_cfg.clearGlobal=1; %last view, clear plotPialSurf from global memory
    else
        sub_cfg.clearGlobal=0;
    end
    sub_cfg_out=plotPialSurf(fsSub,sub_cfg);
    
    % Get electrode colormap limits
    if v==1
        if universalYes(elecCbar)
            elecUsedLimits=sub_cfg_out.elecCbarLimits;
            elecCmapName=sub_cfg_out.elecCmapName;
        end
    end
    if v==1
        if universalYes(olayCbar)
            olayUsedLimits=sub_cfg_out.olayCbarLimits;
            olayCmapName=sub_cfg_out.olayCmapName;
        end
    end
end


%% DRAW COLORBAR(S)
if universalYes(elecCbar) && universalYes(olayCbar),    
    % Colorbar for electrodes
    pos=[.88 .1 .01 .8];
    cbarDGplus(pos,elecUsedLimits,elecCmapName,5,elecUnits);
    
    % Colorbar for pial surface overlay (e.g., neuroimaging)
    pos=[.94 .1 .01 .8];
    cbarDGplus(pos,olayUsedLimits,olayCmapName,5,olayUnits);
elseif universalYes(elecCbar),
    % Colorbar for electrodes
    pos=[.90 .1 .03 .8];
    cbarDGplus(pos,elecUsedLimits,elecCmapName,5,elecUnits);
elseif universalYes(olayCbar)
    % Colorbar for pial surface overlay (e.g., neuroimaging)
    pos=[.90 .1 .03 .8];
    cbarDGplus(pos,olayUsedLimits,olayCmapName,5,olayUnits);
end

% Title
if isfield(cfg,'title')
    if ~isempty(cfg.title)
        % Overall Fig Title
        ht=textsc2014(cfg.title,'title');
        set(ht,'fontweight','bold','fontsize',20,'position',[0.5 0.975]);
    end
end
drawnow;

%subfunction verbReport
function verbReport(report,verbtag,VERBLEVEL)
% verbReport() - Outputs messages if they exceed a desired level of importance.
%
% Usage:
%  >> verbReport(report,verbtag,VERBLEVEL)
%
% Inputs:
%   report    = a string that is some message to the user
%   verbtag   = an intger specifiying the importance of report
%   VERBLEVEL = an integer specifiying a threshold of importance
%               for displaying reports. If verbtag is less than VERBLEVEL, the
%               report will be displayed..
%
% Author: Tom Urbach
% Kutaslab

if nargin<3
    tmpVERBLEVEL = 3;
elseif isempty(VERBLEVEL)
    tmpVERBLEVEL = 3;
else
    tmpVERBLEVEL = VERBLEVEL;
end;

if verbtag <= tmpVERBLEVEL
    if ischar(report)
        fprintf('%s\n',report);
    else
        fprintf('%d\n',report);
    end;
end;
