function downloadHCPatlas(ask_permission)
% function downloadHCPatlas
%
% This function is used to download the HCP atlas and place it into
% the fsaverage label directory. The url to download the atlas files
% is https://figshare.com/ndownloader/articles/3498446/versions/2
%
% Author:
% Noah Markowitz
% Human Brain Mapping Labratory
% NorthShore University Hospital
% October 2021
%

if nargin < 1
    ask_permission = 0;
end

if ask_permission
    prompt = input('Attempting to download the HCP atlas. Enter 1 or "y" to allow and 0 or "n" to exit.  ');
    if universalNo(prompt)
        return;
    end
end
    

% General info
fsDir=getFsurfSubDir();
labelFolder=fullfile(fsDir,'fsaverage','label');

% Download and unzip
url2read = 'https://figshare.com/ndownloader/articles/3498446/versions/2';
outdname_zip = fullfile(labelFolder,'HCP_atlas.zip');
outdname = fullfile(labelFolder,'HCP_atlas');
fprintf('--->Downloading zip file\n');
websave(outdname_zip, url2read);
fprintf('--->Unzipping downloaded file\n');
unzip(outdname_zip, outdname);

% Move files
fprintf('--->Moving files into fsaverage label folder\n');
movefile([outdname  filesep '*.annot'], labelFolder);

% Remove downloaded files and directories
rmdir(outdname,'s');
delete(outdname_zip);

return;
