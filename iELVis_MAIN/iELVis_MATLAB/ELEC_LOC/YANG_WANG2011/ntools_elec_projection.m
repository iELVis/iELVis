function [elec_proj,info_cell] = ntools_elec_projection(ini_loc,ini_pos,row,col,sph,...
    subjectpath,scale,radius,range,std_thres,avg_thres)
% function [elec_proj,info_cell] = ntools_elec_projection(ini_loc,ini_pos,row,col,sph,...
%     subjectpath,scale,radius,range,std_thres,avg_thres)
%
% this program project the calculated grid from ntools_elec_interp_grid to
% the patient's brain lgi outer smoothed surface. 
% 
% Usage: 
% elec_proj = ntools_elec_projection(elec,col,sph,subjectpath);
% elec_proj = ntools_elec_projection(elec,col,sph,subjectpath,scale,radius,range,std_thres,avg_thres)
%
% Inputs:
% elec: cell array with electrodes' names in first column, x y z
%           coordinates in the rest 3
% scale: physical size of dim (in mm)
% col: column number of the grid
% sph: the hemisphere on which the grid locates
% subjectpath: patient's Freesurfer reconstruction folder
% radius: supposed distance between each 2 nearest electrodes (mm)
% range: iteration range
% std_thres: standard deviation (of distance to nearest neighbor) threshold, lower numbers pass (mm)
% avg_thres: average distance to nearest neighbor threshold, the difference of average distance
%                   and radius must be within the this threshold (mm)
%
% Here is how the best projection parameter is chosen:
%1) if any iterations produce a value less than avg thresh, the minimum is taken and std_thresh is ignored
%2) if no iteration produces a value less than avg thresh, then if any iterations produce values less than std thresh the iteration with subtresh std and min avg is taken
%3) If not iteration produces values less than avg thresh or less than std thresh, than you?re supposed to manually an iteration
%
% Outputs:
% elec_proj: cell array contains:
%                   {1}: iteration constant
%                   {2}: average distance of each two nearest electrodes in the projected grid
%                   {3}: standard deviation of the average distances
%                   {4}: center location
%                   {5}: projected grid location
% info_cell: contains other results not being chosen automatically

% default parameter
if ~exist('scale','var') || isempty(scale)
    scale = 1;
end

if ~exist('radius','var') || isempty(radius)
    radius = 10;
end

if ~exist('range','var') || isempty(range)
    range=-5:0.3:-1;
end

if ~exist('std_thres','var') || isempty(std_thres)
    std_thres = 1;
end

if ~exist('avg_thres','var') || isempty(avg_thres)
    avg_thres = 0.1;
end

% read the surf file
surf = fs_read_surf(fullfile(subjectpath,'surf',[sph '.pial-outer-smoothed']));
% This varies depending on the version of FreeSurfer
if ~isfield(surf,'coords')
    if isfield(surf,'vertices')
        surf.coords = surf.vertices;
        surf=rmfield(surf,'vertices');
    else
        surf.coords = surf.coord;
        surf=rmfield(surf,'coord');
    end
end
if size(surf.coords,2)~=3
    surf.coords=surf.coords';
end
if ~isfield(surf,'faces')
    surf.faces=surf.tri;
    surf=rmfield(surf,'tri');
end
if size(surf.faces,2)~=3
    surf.faces=surf.faces';
end

%% get the mesial side data
bb = unique(round(surf.coords(:,2:3)),'rows');
cc = [bb(:,1)-min(bb(:,1))+1 bb(:,2)-min(bb(:,2))+1];

bw = zeros(max(cc(:,1)),max(cc(:,2)));
for i=1:length(cc)
    bw(cc(i,1),cc(i,2)) = 1;
end
bw_close = imclose(bw,strel('square',2));
B = bwboundaries(bw_close); 

bound = [B{1}(:,1)-1+min(bb(:,1)),B{1}(:,2)-1+min(bb(:,2))];
k = dsearchn(surf.coords(:,2:3),bound);
surf_bound = surf.coords(k,:);
surf_bound = unique(surf_bound,'rows');
x = griddata(surf_bound(:,2),surf_bound(:,3),surf_bound(:,1),surf.coords(:,2),surf.coords(:,3));

mesial = (abs(x)-abs(surf.coords(:,1)))>0;

%% interp the grid and get the patch of brain
brain = surf.coords;

% mesial = brain(:,1)<2;
brain(mesial,:) = [];

if length(ini_pos)==2
    [temp_loc,temp_pos] = ntools_elec_223(ini_loc,ini_pos,row,col,sph);
    kk = dsearchn(brain,temp_loc);
    elec = ntools_elec_interp_grid(brain(kk,:),temp_pos,row,col);
elseif length(ini_pos)==3
    elec = ntools_elec_interp_grid(ini_loc,ini_pos,row,col);
elseif length(ini_pos)==4
    kk = dsearchn(brain,ini_loc);
    elec = ntools_elec_locpos4(brain(kk,:), ini_pos, row, col);
end

y = elec(:,3);
z = elec(:,4);
ymin = min(y)-5;
ymax = max(y)+5;
zmin = min(z)-5;
zmax = max(z)+5;
a = brain(:,2)>ymin & brain(:,2)<ymax & brain(:,3)>zmin & brain(:,3)<zmax;
elec_area = brain(a,:);

%% center point and mid point

src = round(rand([round(length(elec_area)/5),1])*(length(elec_area)-1))+1;
ref_pt = elec_area(src,:);

A = zeros(length(ref_pt)-1,3);
b = zeros(length(ref_pt)-1,1);
for i=1:length(ref_pt)-1
    A(i,:) = 2*(ref_pt(i,:)-ref_pt(i+1,:));
    b(i) = sum(ref_pt(i,:).^2-ref_pt(i+1,:).^2);
end

ctr = A\b ;
% ctr = pinv(A)*b ;

% center correction method 1:
%         ctr = [ctr(1)-itr;ctr(2);ctr(3)];

% center correction method 2:
y_mid = (ymax+ymin)/2;
z_mid = (zmax+zmin)/2;
mid_ind =dsearchn(elec_area(:,2:3),[y_mid,z_mid]);
mid = elec_area(mid_ind,:);

%% first loop
info_cell = cell(length(range),5);
fprintf('NO.     itr      distance mean(mm)     standard deviation(mm)     center location         \n');
fprintf('------------------------------------------------------------------------------------------\n');
tic;
for itr = 1:length(range)
    
    ctr1 = ctr+(mid'-ctr)*range(itr);

    [standard, avg, elec_interp] = sphere_proj(ctr1,elec,elec_area,surf,col);

    info_cell(itr,1) = {range(itr)};
    info_cell(itr,2) = {avg*scale};
    info_cell(itr,3) = {standard*scale};
    info_cell(itr,4) = {ctr1};
    info_cell(itr,5) = {elec_interp};

    % fprintf('NO.     itr      distance mean(mm)     standard deviation(mm)     center location         \n');
    fprintf('%.0f    %.2f             %.04f                  %.04f              [ %.0f,%.0f,%.0f  ]      \n',...
        itr,info_cell{itr,1},info_cell{itr,2},info_cell{itr,3},info_cell{itr,4}(1),info_cell{itr,4}(2),info_cell{itr,4}(3));

end

%% second loop

dis_dif = [info_cell{:,2}]-radius; %compare avg distance between projected electrode neighbors with true distance
if ~(all(dis_dif>0) || all(dis_dif<0))
    positive = max(dis_dif,0);
    positive(logical(positive==0)) = nan;
    [~, pos_min] = min(positive);
    negative = min(dis_dif,0);
    negative(logical(negative==0)) = nan;
    [~,neg_max] = max(negative);
    range2_min = min(info_cell{neg_max,1},info_cell{pos_min,1});
    range2_max = max(info_cell{neg_max,1},info_cell{pos_min,1});
    range2 = range2_min+0.05:0.05:range2_max-0.05;
    for itr2 = 1+itr:length(range2)+itr
        
        ctr2 = ctr+(mid'-ctr)*range2(itr2-itr);

        [standard, avg, elec_interp] = sphere_proj(ctr2,elec,elec_area,surf,col);

        info_cell(itr2,1) = {range2(itr2-itr)};
        info_cell(itr2,2) = {avg*scale};
        info_cell(itr2,3) = {standard*scale};
        info_cell(itr2,4) = {ctr2};
        info_cell(itr2,5) = {elec_interp};

        % fprintf('NO.     itr      distance mean(mm)     standard deviation(mm)     center location         \n');
        fprintf('%.0f    %.2f             %.04f                  %.04f              [ %.0f,%.0f,%.0f  ]      \n',...
            itr2,info_cell{itr2,1},info_cell{itr2,2},info_cell{itr2,3},info_cell{itr2,4}(1),info_cell{itr2,4}(2),info_cell{itr2,4}(3));
    end
end

%%
figure(1); clf;
subplot(121);
plot(abs([info_cell{:,2}]-radius),'b-o'); hold on;
axLimits=axis;
plot(axLimits(1:2),[1 1]*avg_thres,'k--');
axLimits=axis;
if axLimits(3)>(avg_thres*.95),
    axis([axLimits(1:2) avg_thres*.95 axLimits(4)]);
elseif axLimits(4)<(std_thres*1.05),
    axis([axLimits(1:3) avg_thres*1.05]);
end
legend('Data','Threshold','location','best');
xlabel('Iteration');
ylabel('mm');
title('Abs(Avg. Dist-True Dist) to Closest Neighbor');

subplot(122);
plot([info_cell{:,3}],'r-o'); hold on;
axLimits=axis;
plot(axLimits(1:2),[1 1]*std_thres,'k--');
axLimits=axis;
if axLimits(3)>(std_thres*.95),
    axis([axLimits(1:2) std_thres*.95 axLimits(4)]);
elseif axLimits(4)<(std_thres*1.05),
    axis([axLimits(1:3) std_thres*1.05]);
end
legend('Data','Threshold','location','best');
xlabel('Iteration');
ylabel('mm');
set(gca,'yaxislocation','right');
title('Std. Distance to Closest Neighbor');

%% output judgement
v = find(abs([info_cell{:,2}]-radius)<=avg_thres);
if isempty(v)
    fprintf('no iteration produced locations with nbor distances lower than the current avg threshold...\n');
    vv = find(cell2mat(info_cell(:,3))<=std_thres);
    if isempty(vv)
        fprintf('no iteration produced locations with nbor distances lower than the current std threshold either, please select the best one from above\n');
        no = input('Please input the NO. of the optimal one from above list, type 0 to quit the program : ');
        if no==0
            error('no data match, quiting the program');
        end
        elec_proj = info_cell(no,:);
    else
        avg_dif = abs([info_cell{vv,2}]-radius);
        tt = find(avg_dif == min(avg_dif));
        elec_proj = info_cell(vv(tt(1)),:);
    end
else
    t = find([info_cell{v,3}]==min([info_cell{v,3}]));
    elec_proj = info_cell(v(t(1)),:);
end
    
toc;

%% subfunc sphere_proj
function  [standard, avg, elec_interp] = sphere_proj(ctr,elec,elec_area,surf,col)

vec1 = repmat(ctr',[size(elec,1) 1])-elec(:,2:4);
vec2 = repmat(ctr',[size(elec_area,1) 1])-elec_area;
for j=1:length(vec1)
    k = dot(repmat(vec1(j,:),[length(vec2) 1]),vec2,2)./(repmat(norm(vec1(j,:)),[length(vec2) 1]).*sqrt(sum(abs(vec2).^2,2)));
    idx(j) = find(k==max(k));
end
elec_interp = elec_area(idx,:);

% find the nearest point and calculate the distance
D0 = dsearchn(surf.coords,elec_interp);
D = fastmarch(surf.faces, surf.coords(:,1), surf.coords(:,2), surf.coords(:,3), D0);
p = 1;
for o = 1:length(D0)
    if (o+1<=length(D0) && mod(o,col)~=0)
        dis1 = D(o,o+1);
        distance(p,:) = [o, o+1, dis1];
        q = p;
    end
    if (o+col<=length(D0))
        dis2 = D(o,o+col);
        distance(q+1,:) = [o, o+col, dis2];
        q = q+1;
    end
    p = q+1;
end

avg = mean(distance(:,3));
standard = std(distance(:,3));

