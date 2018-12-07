function repoVersion=iELVis_getGitInfo()
% function repoVersion=iELVis_getGitInfo()
%
% Returns the version of the iELVis repo being used (i.e., the git hash)
% This function assumes that repo root directory is wherever
% iELVis_getGitInfo.m is.

fullfname=which('iELVis_getGitInfo.m');
[filepath,name,ext] = fileparts(fullfname);
repoVersion='';
% Just in case of weirdness across OSes
try
    gitInfo=getGitInfo(filepath);
    repoVersion=gitInfo.hash;
catch
    
end
