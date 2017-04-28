function elec = ntools_elec_calc_depth(ini_cell)

% calculate the depth electrodes with initial locations, ini_cell is a cell
% array in which first column is the [elec_name number], the rest 3
% columns are the x y z coordinates
%
% output: elec is a structure contains the depth names and their locations,
% data are only the coordinates without depth names, used for binary image
% conversion

if isempty(ini_cell)
    disp('No depth electrode.');
    elec = [];
    return;
end
fprintf('Calculating the depth electrodes....'); tic;
name = regexp(ini_cell(:,1),'[A-Za-z]*[^\d*]','match');
for i=1:length(name)
    ini_name(i) = name{i};
end
name = unique(ini_name);

for i = 1:length(name)
        elec_temp = cell(size(ini_cell));
%         n = strfind(ini_depth(:,1),name{i});
        n = regexp(ini_cell(:,1),[name{i} '[\d]'],'match');
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
        if length(elec_num)~=2
            error('only 2 initial positions are required');
        end
        elec_ini_loc = elec_temp(:,2:4);
        [E_temp, E_temp_pos]= ntools_elec_locpos(cell2mat(elec_ini_loc),[str2double(cell2mat(elec_num{1}));str2double(cell2mat(elec_num{2}))]);
        elec_temp2.(char(name{i})) = [E_temp_pos E_temp];

end

%% get the data from the struct
l = 1;
for i = 1:length(name)
    for j = 1:size(elec_temp2.(char(name{i})),1)
        elec_name(l) = cellstr(sprintf('%s%.2d',char(name{i}),elec_temp2.(char(name{i}))(j,1)));
        elec_pos(l,:) = elec_temp2.(char(name{i}))(j,2:4);
        l = l+1;
    end
end
elec = cell(size(elec_name,2),5);

elec(:,1) = upper(elec_name)';
elec(:,2:4) = num2cell(elec_pos);

elec(:,5) = repmat({'D'},[l-1 1]);

fprintf('Done. (%f seconds) \n\n', toc);