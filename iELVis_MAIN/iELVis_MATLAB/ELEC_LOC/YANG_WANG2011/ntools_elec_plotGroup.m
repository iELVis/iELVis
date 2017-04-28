function ntools_elec_plotGroup(varargin)

% ntools_elec_plot(elec_file,surf_file,'labelshow',0/1,'genimg',0/1)
%
% a stand-alone program that shows ieeg electrodes on the pial surface and
% save the images into textfile folder/images/. Default saving name is
% PatientID(NYxx)_space_elecLabel_viewpoint_hemisphere.png
% 
% space could be T1 or MNI
% elecLabel: numbers in the input elec_file, e.g. 1,2,3,4,5. Type the
%           group name when the script asks for description
% viewpoint: could be left, right, top, below, back or front
% hemisphere: could be lh, rh, or both
% 
% required input:
% elec_text: text or xls(x) file with xyz electrode coordinates
% pial_mat: matlab structure with pial surface
% plot: elecLabel to plot, e.g. 2,4
% labelshow: to show the electrode labels (1) or not (0)
% genimg: to save the images (1) or not (0)
% groupname: How do you want to describe these electrodes, e.g.
%           'seizure onset', or leave it blank ([])
%
% Usage: run ntools_elec_plot in command window
% the gui is prompted out for file selection
% or: ntools_elec_plot(textfilename, {'lh.mat', 'rh.mat'})
%     ntools_elec_plot(textfilename,lh.mat,'plot',[2,3,9],'labelshow',1,'genimg',1,'groupname','seizure onset');
%
% by  Hugh Wang, Xiuyuan.Wang@nyumc.org, Aug 13th, 2013

%% Get the elec info
if nargin==0
    [FileName,PathName] = uigetfile({'*.xlsx';'*.txt';'*.xls'},'Select the electrodes text file',pwd); 
    [surfname, surfpath] = uigetfile('*.mat','Select the patient brain surf',PathName,'MultiSelect','on');
    surf = strcat(surfpath,surfname);      
elseif nargin>=2
    aa = strfind(varargin{1},'/');
    if isempty(aa)
        PathName = pwd;
        FileName = varargin{1};
    else
        FileName = varargin{1}(aa(end)+1:end);
        PathName = varargin{1}(1:aa(end));
    end
    surf = varargin{2}; 

    try labelshow = varargin{find(strcmp('labelshow',varargin))+1}; catch err; end
    try genimg = varargin{find(strcmp('genimg',varargin))+1}; catch err; end
    try showpart = varargin{find(strcmp('groupname',varargin))+1}; catch err; end
    try plt = varargin{find(strcmp('plot',varargin))+1}; catch err; end
end

if exist(fullfile(PathName, FileName),'file')
    [~,~,ext] = fileparts(FileName);
    if strcmpi(ext,'.txt')
        fid = fopen(fullfile(PathName, FileName));
        elec_all = textscan(fid,'%s %f %f %f %f');
        elec_cell = [elec_all{1},num2cell(elec_all{2}),num2cell(elec_all{3}),num2cell(elec_all{4})];
        fclose(fid);
    elseif strcmpi(ext,'.xls') || strcmpi(ext,'.xlsx')
        [~,~,elec_all] = xlsread(fullfile(PathName,FileName)); 
        elec_cell = [elec_all{1},num2cell(elec_all{2}),num2cell(elec_all{3}),num2cell(elec_all{4})];
    end
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

%% Separate the elecs by category numbers
catnum = elec_all{5};
uninum = unique(catnum,'sorted');
cmap = jet(length(uninum));


%% Plot the elecs
if ~exist('plt','var')
    plt = input(sprintf('What category do you want to plot? (%s~%s): ',num2str(uninum(1)),num2str(uninum(end))),'s');
    plt = str2num(plt);
end

if ~exist('showpart','var')
    showpart = input('How do you want to describe these electrodes? (e.g. seizure onset): ','s');
end
showpart = regexprep(showpart,' ','_');

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

nyu_plot(surf_brain,sph,[],[]);

for i=1:length(plt)
    idx = find(catnum==plt(i));
    elec = cell2mat(elec_cell(idx,2:4));
    elec_name = char(elec_cell(idx,1));
    col = cmap(uninum==plt(i),:);
    for j=1:size(elec,1)
        plotSpheres(elec(j,1),elec(j,2),elec(j,3),2,col);
        if labelshow==1
            [xx, yy, zz] = adjust_elec_label(elec(j,:)); % default radius = 2
            text('Position',[xx yy zz],'String',elec_name(j,:),'Color','w','VerticalAlignment','top');
        end
    end 
    hold on;
end

colormap(cmap);
cbar = colorbar('YTickLabel',{uninum},'Clim',[uninum(1),uninum(end)]);
ytick_range = get(cbar,'YLim');
interval = (ytick_range(2)-ytick_range(1))/(length(uninum));
ytick = [ytick_range(1):interval:ytick_range(2)]+interval/2;
set(cbar,'YTick',ytick);

hold off;

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

%%
function [shand]=plotSpheres(spheresX, spheresY, spheresZ, spheresRadius,varargin)

if nargin>4,
    col=varargin{:};
end

spheresRadius = ones(length(spheresX),1).*spheresRadius;
% set up unit sphere information
numSphereFaces = 25;
[unitSphereX, unitSphereY, unitSphereZ] = sphere(numSphereFaces);

sphereCount = length(spheresRadius);

% for each given sphere, shift the scaled unit sphere by the
% location of the sphere and plot
for i=1:sphereCount
sphereX = spheresX(i) + unitSphereX*spheresRadius(i);
sphereY = spheresY(i) + unitSphereY*spheresRadius(i);
sphereZ = spheresZ(i) + unitSphereZ*spheresRadius(i);
%shand=surface(sphereX, sphereY, sphereZ,'FaceColor',col,'EdgeColor',col,'Tag','btq_sphere');
shand=surface(sphereX, sphereY, sphereZ,'FaceColor',col,'EdgeColor','none','AmbientStrength',0.7);
end
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