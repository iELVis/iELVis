function ntools_elec_plot(varargin)

% a stand-alone program that shows ieeg electrodes on the pial surface and
% save the images into textfile folder/images/. Default saving name is
% PatientID(NYxx)_space_elecLabel_viewpoint_hemisphere.png
% 
% space could be T1 or MNI
% elecLabel could be grid, strip, depth, or grid and strip
% viewpoint could be left, right, top, below, back or front
% hemisphere could be lh, rh, or both
% 
% required input:
% elec_text: text file with xyz electrode coordinates
% pial_mat: matlab structure with pial surface
%
% optional input:
% plt: part to plot 
%     1 Grid only 
%     2 Strip only 
%     3 Depth Only 
%     4 Both Grid and Strip  
% labelshow: to show the electrodes' labels (1) or not (0)
% genimg: to save the images (1) or not (0)
% parcellation: to save the image with freesurfer parcellation (1) or
% not (0)
%
% Usage: run ntools_elec_plot in command window
% the gui is prompted out for file selection
% or: ntools_elec_plot(textfilename, {'lh.mat', 'rh.mat'})
%       ntools_elec_plot(textfilename, l(r)h.mat);
%
% also see: http://ieeg.pbworks.com/Viewing-Electrode-Locations
%
% written by  Hugh Wang, Xiuyuan.Wang@nyumc.org, May 13, 2009

% modified on May 14, 2009 by Hugh
% make judgement on the input file type and not sensitive to the order of 
% input variable. add the option to show the electrodes' labels or not.

% modified on July 22nd, 2009 by Hugh
% for subjects who has electrodes on both hemisphere, loading the both
% pial.mat file will generate the image with whole brain and save the
% images from 6 views (left, right, top, below, back & front)

% modified on Aug 8th, 2009 by Hugh
% show only lh(rh)'s electrodes if choose lh(rh)_surf.mat

% modified on Jan 28th, 2010 by Hugh
% default saving name

%% Get the elec info
if nargin==0
    [FileName,PathName] = uigetfile('*.txt','Select the electrodes text file','/home/halgdev/projects/nyuproj/loc/'); % NYU settings
    [surfname, surfpath] = uigetfile('*.mat','Select the patient brain surf',PathName,'MultiSelect','on');
    surf = strcat(surfpath,surfname);      
elseif nargin>=2
    aa = strfind(varargin{1},'/');
    FileName = varargin{1}(aa(end)+1:end);
    PathName = varargin{1}(1:aa(end));
    surf = varargin{2}; 

    try labelshow = varargin{find(strcmp('labelshow',varargin))+1}; catch err; end
    try genimg = varargin{find(strcmp('genimg',varargin))+1}; catch err; end

end

if exist([PathName, FileName],'file')
    fid = fopen([PathName, FileName]);
    elec_all = textscan(fid,'%s %f %f %f %s');
    elec_cell = [elec_all{1},num2cell(elec_all{2}),num2cell(elec_all{3}),num2cell(elec_all{4})];
else
    elec_cell = [];
end

%% Get the filename info
b = strfind(FileName,'_');
Pname = FileName(1:b(1)-1);

if ~isempty(strfind(upper(FileName),'T1'))
    space = '_T1_';
elseif ~isempty(strfind(upper(FileName),'MNI'))
    space = '_MNI_';
else
    space = '_';
end

if length(surf)==2
    sph = 'both';
else
    sph = regexp(lower(surf),'[r,l]h','match');
    sph = char(sph{:});
end

%% Separate Grid, Strip and Depth electrodes

if isempty(char(elec_all{5}(:)))
    g = strmatch('G',upper(elec_cell(:,1)));
    d = strmatch('D',upper(elec_cell(:,1)));
else
    g = strmatch('G',upper(elec_all{5}));
    d = strmatch('D',upper(elec_all{5}));
end

if ~isempty(g) && ~isempty(d)
    elec_grid = elec_cell(g,:);
    elec_depth = elec_cell(d,:);
    elec_cell([g;d],:) = [];    
elseif isempty(d)
    elec_depth = [];
    elec_grid = elec_cell(g,:);
    elec_cell(g,:) = [];
elseif isempty(g)
    elec_grid = [];
    elec_depth = elec_cell(d,:);
    elec_cell(d,:) = [];
end


%% Plot the elecs
if ~exist('plt','var')
    plt = menu('What part do you want to plot?','Grid only', 'Strip only','Depth Only','Both Grid and Strip','Brain only');
end
if ~exist('labelshow','var')
    labelshow = menu('Do you want to show the label?','Yes','No');
end
if ~exist('genimg','var')
    genimg = menu('Do you want to save the images?','Yes', 'No');
end

if strcmp(sph,'both')
    surf_brain.sph1 = load(surf{1});
    surf_brain.sph2 = load(surf{2});
else 
    surf_brain = load(surf);
end

if plt==1 && ~isempty(elec_grid)
    showpart = 'G';
    nyu_plot(surf_brain,sph,cell2mat(elec_grid(:,2:4)),char(elec_grid(:,1)),'r',labelshow);
elseif plt==2 && ~isempty(elec_cell)
    showpart = 'S';
    nyu_plot(surf_brain,sph,cell2mat(elec_cell(:,2:4)),char(elec_cell(:,1)),'b',labelshow);
elseif plt==3 && ~isempty(elec_depth)
    showpart = 'D';
    nyu_plot(surf_brain,sph,cell2mat(elec_depth(:,2:4)),char(elec_depth(:,1)),'g',labelshow,1.5,0.3);
elseif plt==4 && ~isempty(elec_grid) && ~isempty(elec_cell)
    showpart = 'GS';
    elec = cell2mat(elec_cell(:,2:4));
    elec_name = char(elec_cell(:,1));
    nyu_plot(surf_brain,sph,cell2mat(elec_grid(:,2:4)),char(elec_grid(:,1)),'r',labelshow); hold on;
    for i=1:size(elec,1)
        plotSpheres(elec(i,1),elec(i,2),elec(i,3),2,'b');
        if labelshow==1
            [xx, yy, zz] = adjust_elec_label(elec(i,:)); % default radius = 2
            text('Position',[xx yy zz],'String',elec_name(i,:),'Color','w','VerticalAlignment','top');
        end
    end
    hold off;  
elseif plt==5
    showpart = 'Brain';
    nyu_plot(surf_brain,[],[],[]);
else
    disp('sorry, the electrodes you choose to show are not on the surface you loaded');
    return;
end

%% save images

if genimg==1
    if ~exist([PathName 'images/'],'dir')
        mkdir([PathName 'images/']);
    end
    
    if labelshow==1
        label = '_label';
    else
        label = [];
    end
    
    format = 'png';
    if strcmp(sph,'lh')
        view(270, 0);
        saveas(gcf,[PathName,'images/',Pname,space,showpart,'_lateral_',sph,label],format);
        view(90,0);
        saveas(gcf,[PathName,'images/',Pname,space,showpart,'_mesial_',sph,label],format);
        
    elseif strcmp(sph,'rh')
        view(270, 0);
        saveas(gcf,[PathName,'images/',Pname,space,showpart,'_mesial_',sph,label],format);
        view(90,0);
        saveas(gcf,[PathName,'images/',Pname,space,showpart,'_lateral_',sph,label],format);
        
    elseif strcmp(sph,'both')
        view(270, 0);
        saveas(gcf,[PathName,'images/',Pname,space,showpart,'_left_',sph,label],format);
        view(90,0);
        saveas(gcf,[PathName,'images/',Pname,space,showpart,'_right_',sph,label],format);
    end
    view(0,0);
    saveas(gcf,[PathName,'images/',Pname,space,showpart,'_posterior_',sph,label],format);

    view(180,0);
    saveas(gcf,[PathName,'images/',Pname,space,showpart,'_frontal_',sph,label],format);

    view(90,90);
    saveas(gcf,[PathName,'images/',Pname,space,showpart,'_dorsal_',sph,label],format);

    view(90,-90);
    set(light,'Position',[1 0 -1]);
    saveas(gcf,[PathName,'images/',Pname,space,showpart,'_ventral_',sph,label],format);
else 
    return;
end

end

%% subfunction 
%% nyu_plot
function nyu_plot(surf_brain,sph,elec,elecname,color,label,radius,alpha)

if ~exist('color','var')
    color = 'w';
end
if ~exist('label','var')
    label = 0;
end
if ~exist('alpha','var')
    alpha = 1;
end
if ~exist('radius','var')
    radius = 2;
end

figure;

col = [.7 .7 .7];
if strcmp(sph,'both')
    sub_sph1.vert = surf_brain.sph1.coords;
    sub_sph1.tri = surf_brain.sph1.faces;

    sub_sph2.vert = surf_brain.sph2.coords;
    sub_sph2.tri = surf_brain.sph2.faces;
    
    col1=repmat(col(:)', [size(sub_sph1.vert, 1) 1]);
    col2=repmat(col(:)', [size(sub_sph2.vert, 1) 1]);
    
    trisurf(sub_sph1.tri, sub_sph1.vert(:, 1), sub_sph1.vert(:, 2),sub_sph1.vert(:, 3),...
        'FaceVertexCData', col1,'FaceColor', 'interp','FaceAlpha',alpha);
    hold on;
    trisurf(sub_sph2.tri, sub_sph2.vert(:, 1), sub_sph2.vert(:, 2), sub_sph2.vert(:, 3),...
        'FaceVertexCData', col2,'FaceColor', 'interp','FaceAlpha',alpha);
else    
    if isfield(surf_brain,'coords')==0
        sub.vert = surf_brain.surf_brain.coords;
        sub.tri = surf_brain.surf_brain.faces;
    else
        sub.vert = surf_brain.coords;
        sub.tri = surf_brain.faces;
    end
    col=repmat(col(:)', [size(sub.vert, 1) 1]);
    trisurf(sub.tri, sub.vert(:, 1), sub.vert(:, 2), sub.vert(:, 3),...
        'FaceVertexCData', col,'FaceColor', 'interp','FaceAlpha',alpha);
end

shading interp;
lighting gouraud;
material dull;
light;
axis off
hold on;
for i=1:size(elec,1)
    plotSpheres(elec(i,1),elec(i,2),elec(i,3),radius,color);
    if label==1
        [x, y, z] = adjust_elec_label(elec(i,:),radius);
        text('Position',[x y z],'String',elecname(i,:),'Color','w','VerticalAlignment','top');
    end
end
set(light,'Position',[-1 0 1]); 
    if strcmp(sph,'lh')
        view(270, 0);      
    elseif strcmp(sph,'rh')
        view(90,0);        
    elseif strcmp(sph,'both')
        view(90,90);
    end
set(gcf, 'color','black','InvertHardCopy', 'off');
axis tight;
axis equal;
end

%% adjust_elec_label
function [x, y, z] = adjust_elec_label(elec,radius)

if ~exist('radius','var')
    radius = 2;
end

if elec(1)>0
    x = elec(1)+radius;
else
    x = elec(1)-radius;
end

if elec(3)>0
    z = elec(3)+radius;
else
    z = elec(3)-radius;
end

y = elec(2);

end

%% plotSpheres
function [shand]=plotSpheres(spheresX, spheresY, spheresZ, spheresRadius,varargin)

if nargin>4,
    col=varargin{:};
end

spheresRadius = ones(length(spheresX),1).*spheresRadius;
% set up unit sphere information
numSphereFaces = 25;
[unitSphereX, unitSphereY, unitSphereZ] = sphere(numSphereFaces);

% set up basic plot
sphereCount = length(spheresRadius);

% for each given sphere, shift the scaled unit sphere by the
% location of the sphere and plot
for i=1:sphereCount
sphereX = spheresX(i) + unitSphereX*spheresRadius(i);
sphereY = spheresY(i) + unitSphereY*spheresRadius(i);
sphereZ = spheresZ(i) + unitSphereZ*spheresRadius(i);
shand=surface(sphereX, sphereY, sphereZ,'FaceColor',col,'EdgeColor','none','AmbientStrength',0.7);
end

end
