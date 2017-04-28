function [elec, elec_stats, info_cell, data]= ntools_elec_calc_grid(ini_cell, subjectpath,scale,radius)

% calculate the grid electrodes with initial locations, ini_cell is a cell
% array in which first column is the [elec_name number], the rest 3
% columns are the x y z coordinates
%
% output: elec is a structure contains the grid names and their locations,
% data are only the coordinates without grid names, used for binary image
% conversion

if isempty(ini_cell)
    disp('No grid electrode.');
    elec = []; data = []; elec_stats = []; info_cell = [];
    return;
end

fprintf('Calculating the grids.....\n')

name = regexp(ini_cell(:,1),'[A-Za-z]*[^\d*]','match');
for i=1:length(name)
    ini_name(i) = name{i};
end
name = unique(ini_name); %# of unique grids

% get the grid initial positions by name
for i = 1:length(name)
    elec_temp = cell(size(ini_cell));
    n = regexp(ini_cell(:,1),['^' name{i} '[\d]'],'match'); % # of electrodes in this grid
    k = 1;
    for l = 1:length(n)
        if ~isempty(n{l})
            elec_temp(k,:) = ini_cell(l,:);
            k = k+1;
        end
    end
    elec_temp(all(cellfun(@isempty,elec_temp),2),:) = [];
    elec_num = regexp(elec_temp(:,1),'[^A-Za-z]*[\d*]','match');
    elec_num(all(cellfun(@isempty,elec_num),2),:) = [];
    
    for ll = 1:length(elec_num)
        ini_pos(ll) = str2double(cell2mat(elec_num{ll})); % electrode #s
    end
    ini_loc = cell2mat(elec_temp(:,2:4)); % electrode initial RAS coordinates

    % determine the hemisphere that grid locates
    if ini_loc(:,1)>0
        sph = 'rh';
    elseif ini_loc(:,1)<0
        sph = 'lh';
    else
        error(['wrong initial positions for grid ', name{i}]);
    end
    
    % input the size of the grid
    a = menu(['What is the size of the grid ',name{i},' ?'], '8*8', '4*8', '4*5','2*8','manually input');
    if a==1
        s = [8,8];
    elseif a==2
        s = [4,8];
    elseif a==3
        s = [4,5];
    elseif a==4
        s = [2,8];
    else
        s = input(['Please input the size of the grid ',name{i},' [row column]: ']);
    end
    
    % project the grid on the outer-brain surface
    [elec_proj, info_cell]= ntools_elec_projection(ini_loc,ini_pos,s(1),s(2),sph,subjectpath,scale,radius);
    elec_temp2.(char(name{i})) = elec_proj{1,5};
    elec_stats(i,:) = elec_proj(1,1:4);
    
    % clear the unnecessary data
    clear elec_temp elec_num elec_proj ini_pos ini_loc;
	fprintf('Grid stats: itr#: %f mean: %f std:%f\n\n',cell2mat(elec_stats(i,1:3)))
end

%% get the data from the struct
l = 1;
for i = 1:length(name)
    for j = 1:size(elec_temp2.(char(name{i})),1)
        name_num(l) = cellstr(sprintf('%s%.2d',char(name{i}),j));
        elec_pos(l,:) = elec_temp2.(char(name{i}))(j,:);
        l = l+1;
    end
end

elec = cell(size(name_num,2),5);
elec(:,1) = upper(name_num)';
elec(:,2:4) = num2cell(elec_pos);

elec(:,5) = repmat({'G'},[l-1 1]);

fprintf('Done \n\n');