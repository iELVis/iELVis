%%
eNamesFile='/Applications/freesurfer/subjects/TWH011/elec_recon/TWH011_left.electrodeNames';
eCoordFile='/Applications/freesurfer/subjects/TWH011/elec_recon/TWH011_left.PIALVOX';
eNames=csv2cell(eNamesFile);
eCoords=dlmread(eCoordFile);
ini_cell=[];
ct=0;
for a=1:length(eNames),
    if ~isempty(findstr('Grid',eNames{a}))
        ct=ct+1;
        
        splitId=find(eNames{a}=='_');
        ini_cell{ct,1}=[eNames{a}(1:(splitId-1)) ' ' eNames{a}((splitId+1):end)];
        %ini_cell{ct,2:4}=eCoords(a,2:4);
        ini_cell{ct,2}=eCoords(a,2);
        ini_cell{ct,3}=eCoords(a,3);
        ini_cell{ct,4}=eCoords(a,4);
    end
end



%%
subjectpath='/Applications/freesurfer/subjects/TWH011';
scale=1; % not sure what this should be
%scale: physical size of dim (in mm)<-default is 1, which I believe means
%that you are feeding it coordinates in units of mm

radius=10; % used by ntools_elec_projection(ini_loc,ini_pos,s(1),s(2),sph,subjectpath,scale,radius);
% radius: supposed distance between each 2 nearest electrodes (mm)

%%
%function [elec, elec_stats, info_cell, data]= ntools_elec_calc_grid(ini_cell, subjectpath,scale,radius)

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
name = unique(ini_name);

% get the grid initial positions by name
for i = 1:length(name) % length(name)=# of different grid electrodes
    elec_temp = cell(size(ini_cell));
    n = regexp(ini_cell(:,1),['^' name{i} '[\d]'],'match');
    k = 1;
    for l = 1:length(n)
        if ~isempty(n{l})
            elec_temp(k,:) = ini_cell(l,:);
            k = k+1;
        end
    end
    elec_temp(all(cellfun(@isempty,elec_temp),2),:) = []; % subset of ini_cell corresponding to this grid
    elec_num = regexp(elec_temp(:,1),'[^A-Za-z]*[\d*]','match'); % The numerical index of each grid electrod
    elec_num(all(cellfun(@isempty,elec_num),2),:) = [];
    
    for ll = 1:length(elec_num)
        ini_pos(ll) = str2double(cell2mat(elec_num{ll}));
    end
    ini_loc = cell2mat(elec_temp(:,2:4));

    % determine the hemisphere that grid locates
    if ini_loc(:,1)>0
        sph = 'rh';
    elseif ini_loc(:,1)<0
        sph = 'lh';
    else
        error(['wrong initial positions for grid ', name{i}]);
    end
    sph='lh'; % ?? DG fix later
    
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