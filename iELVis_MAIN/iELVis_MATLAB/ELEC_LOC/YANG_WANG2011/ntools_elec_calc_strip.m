function [elec, data] = ntools_elec_calc_strip(ini_cell,subjectpath,sph)
%function [elec, data] = ntools_elec_calc_strip(ini_cell,subjectpath,sph)
%
% find the strip electrodes on the smoothed pial surface using nearest points

% History:
% 2015-10: updated to deal with columns added to ini_cell (chan # & hem)

if isempty(ini_cell)
    disp('no strip electrodes detected');
    elec = []; data = [];
    return;
end
fprintf('Calculating the strip electrodes....'); tic;
%strip = cell2mat(ini_cell(:,2:4));
strip = cell2mat(ini_cell(:,3:5));

if strcmp(sph,'both')
    surf_lh = fs_read_surf([subjectpath '/surf/lh.pial-outer-smoothed']);
    if ~isfield(surf_lh,'coords')
        surf_lh.coords = surf_lh.vertices;
    end
    surf_rh = fs_read_surf([subjectpath '/surf/rh.pial-outer-smoothed']);
    if ~isfield(surf_rh,'coords')
        surf_rh.coords = surf_rh.vertices;
    end    
    surf = [surf_lh.coords;surf_rh.coords];
else
    surf_h = fs_read_surf([subjectpath '/surf/',sph,'.pial-outer-smoothed']);
    if ~isfield(surf_h,'coords')
        if isfield(surf_h,'coord')
            % Apparentlly this varies with freesurfer version
            surf_h.coords=surf_h.coord;
        else
            surf_h.coords = surf_h.vertices;
        end
    end
    surf = surf_h.coords;
end
if size(surf,1)==3,
    surf=surf'; % Apparently this varies with freesurfer version
end
k = dsearchn(surf,strip);
data = surf(k,:);

% data = ICP_finite(surf,strip,struct('Optimizer','fminsearch'));
%ini_cell(:,1) = regexprep(ini_cell(:,1),'(?<!\d)(\d)(?!\d)','0$1');
%elec = [ini_cell(:,1), num2cell(data), repmat({'S'},[length(strip),1])];
elec = [ini_cell(:,1:2), num2cell(data), ini_cell(:,6:7)];

fprintf('Done. (%f seconds) \n\n', toc);