function coord = depths_interpol(fs_subj)
%
% interpolate coordinates of stereotactic electrodes based on contacts
% located at the extremities of the shaft.
% 
% Please note that this function belong to the iElvis toolboxe
% and is therefore subjected to the same regulations
% (contacts: david.m.groppe@gmail.com ; manuel.mercier@a3.epfl.ch)
%
% input: fs_subj name of the subject in the FS folder
%
% output: coord structure containing 
%           >> elec: {nx1 cell}
%           >> location: {nx1 cell}
%
% with:
% - elec: electrodes name
% - location: brain region where the centroid of each electrode is localized based on freesurfer parcellation
%
%
% files needed:
% - electrodes coordinates (*.3dUndump.VOX)
% - electrodes name (*.electrodeNames)
% 
% 09-2017: created by MrM;
%

fsDir=getFsurfSubDir();

% load electrodes names & Voxel coordinates
coordFname=fullfile(fsDir,fs_subj,'elec_recon',sprintf('%sPostimpLoc.txt',fs_subj));
copyfile(coordFname, [label_file(1:end-4) '_beforeDepthsInterpolation.txt']);

fid=fopen(coordFname);
IO_mgrid=textscan(fid,'%s %s %s %s %s %s %s');
fclose(fid);
for i=1:length(IO_mgrid{1})
    label{i,1} = strcat(IO_mgrid{1}{i},IO_mgrid{2}{i},'_',IO_mgrid{6}{i},'_',IO_mgrid{7}{i});
    coord(i,:)=[str2num(IO_mgrid{3}{i}),str2num(IO_mgrid{4}{i}),str2num(IO_mgrid{5}{i})];
end


%% get extremities of electrode shaft

% extract electrode shaft names
e_shaft = regexp(label,'\D+(\.)?(\D+)?','match');
for i=1:length(e_shaft)
   shaft_name{i} = strcat(e_shaft{i}{1,1},e_shaft{i}{1,2});
end
clear e_shaft;
[shaft_name, idx] = unique(shaft_name, 'stable');
shaft_info = cell2struct(shaft_name,'name');
clear shaft_name

idx(1:end-1,2) = idx(2:end,1)-1;
idx(end,2) = length(label);
shaft_ext_1 = coord(idx(:,1),:);
shaft_ext_2 = coord(idx(:,2),:);

%% get electrodes shaft parameters

for i=1:length(shaft_info)
    shaft_info(i).interpol = input(['Do you want to interpolate ' shaft_info(i).name ' ? (yes=1;no=0)\n' ]);
    if universalYes(shaft_info(i).interpol)
        shaft_info(i).nbElect   = input(['how many electrodes are on the shaft ' shaft_info(i).name ' ? \n' ]);
    else
        shaft_info(i).nbElect = NaN;
    end
end

%% Interpolate
n=0;
for i=1:length(shaft_info)
    if shaft_info(i).interpol
        n=n+1;
        for xyz=1:3
            IO_mgrid{2+xyz}(idx(i,1):idx(i,2)) = num2cell(linspace(shaft_ext_1(i,xyz),shaft_ext_2(i,xyz),shaft_info(i).nbElect)');
        end
    end
end

%% Save

save(fullfile(fsDir,fs_subj,'elec_recon\Interpolation_Info.mat'),'shaft_info');

fid=fopen(coordFname,'w');
for i=1:size(IO_mgrid{1},1)
    fprintf(fid,'%s %s %s %s %s %s %s\n', ....
        IO_mgrid{1}{i},IO_mgrid{2}{i},IO_mgrid{3}{i},IO_mgrid{4}{i},IO_mgrid{5}{i},IO_mgrid{6}{i},IO_mgrid{7}{i});
end
fclose(fid);


end