function depths_mergeshafts(fs_subj)

% In BioImage Suite, the number of contacts per electrode shaft is limited to 16.
% Longuer shaft can be managed as two smaller ones,
% corresponding names should contains a suffix to indicate their location from the opening in the skull
% (the furthest is #1).
%
% This function merge these shaft(s) by:
% - renaming the shaft without the suffix number,
% - renumbering the contacts,
% - reorganizing the order of the electrode accordingly to the input file 
%
% It uses the output of the function "makeIniLocTxtFile" from iElvis toolbox
% 
% 08-2017: created by MrM;
%


fs_dir=getFsurfSubDir();
% Folder with surface files
recon_folder=fullfile(fs_dir,fs_subj,'elec_recon\');

% load electrodes name and make a copy
files=dir([recon_folder fs_subj 'PostimpLoc.txt']);
n=size(files,1);
if n==1
    label_file=fullfile(recon_folder,files.name);
elseif n==0
    disp('No electrodeNames file found. Please do it manualy.');
    [temp_file,elec_dir]=uigetfile([recon_folder '*.electrodeNames'],'select electrode names file');
    label_file=fullfile(elec_dir,temp_file);
    clear elec_dir temp_file files n
elseif n>1
    disp('More than one electrodeNames file found. Please do it manualy.');
    [temp_file,elec_dir]=uigetfile([recon_folder '*.electrodeNames'],'select electrode names file');
    label_file=fullfile(elec_dir,temp_file);
    clear elec_dir temp_file files n
end
copyfile(label_file, [label_file(1:end-4) '_Copy.txt'])

fid=fopen(label_file);
import=textscan(fid,'%s %d %f %f %f %s %s\n');
fclose(fid);
suffix = regexp(import{1},'\d+(\.)?(\d+)?','match');
idx = find(~cellfun(@isempty,suffix) & [0; ~cellfun(@strcmp,import{1}(1:end-1),import{1}(2:end))]);

export = import;
for i=1:2:length(idx)
    export{1}(idx(i):idx(i+1)+import{2}(idx(i+1))-1) = regexp(import{1}{idx(1)},'\D+(\.)?(\D+)?','match');
    if str2num(cell2mat(suffix{idx(i)})) < str2num(cell2mat(suffix{idx(i+1)}))
        part1 = idx(i):(idx(i)+import{2}(idx(i))-1);
        part2 = idx(i+1):(idx(i+1)+import{2}(idx(i+1))-1);
        if import{2}(idx(i)) < import{2}(idx(i)+1) && import{2}(idx(i+1)) < import{2}(idx(i+1)+1) 
           export{2}([part1,part2]) = 1:(import{2}(idx(i))+ import{2}(idx(i+1)));
        elseif import{2}(idx(i)) > import{2}(idx(i)+1) && import{2}(idx(i+1)) > import{2}(idx(i+1)+1)   
            export{2}([part1,part2]) = (import{2}(idx(i))+ import{2}(idx(i+1))):-1:1;
            export{3}(part2) = import{3}(part1);
            export{4}(part2) = import{4}(part1);
            export{5}(part2) = import{5}(part1);
            export{3}(part1) = import{3}(part2);
            export{4}(part1) = import{4}(part2);
            export{5}(part1) = import{5}(part2);
        end
            
    elseif str2num(cell2mat(suffix{idx(i)})) > str2num(cell2mat(suffix{idx(i+1)}))
        part1 = idx(i+1):(idx(i+1)+import{2}(idx(i+1))-1);            
        part2 = idx(i):(idx(i)+import{2}(idx(i))-1);
        if import{2}(idx(i)) < import{2}(idx(i)+1) && import{2}(idx(i+1)) < import{2}(idx(i+1)+1)            
            export{2}([part2,part1]) = 1:(import{2}(idx(i))+ import{2}(idx(i+1)));
            export{3}(part2) = import{3}(part1);
            export{4}(part2) = import{4}(part1);
            export{5}(part2) = import{5}(part1);
            export{3}(part1) = import{3}(part2);
            export{4}(part1) = import{4}(part2);
            export{5}(part1) = import{5}(part2);
        elseif import{2}(idx(i)) > import{2}(idx(i)+1) && import{2}(idx(i+1)) > import{2}(idx(i+1)+1)
           export{2}([part2,part1]) = (import{2}(idx(i))+ import{2}(idx(i+1))):-1:1;
        end
    end
end

fid=fopen(label_file,'w');
for i=1:length(export{1})
    fprintf(fid,'%s %d %f %f %f %s %s\n', ...
        export{1}{i},export{2}(i),export{3}(i),export{4}(i),export{5}(i),export{6}{i},export{7}{i});
end
fclose(fid);

end