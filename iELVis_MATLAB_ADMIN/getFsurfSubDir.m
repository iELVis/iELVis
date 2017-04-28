function fsurfSubDir = getFsurfSubDir()
%function fsurfSubDir = getFsurfSubDir()
%   Returns the Freesurfer Subject Directory as defined by the global
%   Matlab variable globalFsDir or the shell variable SUBJECTS_DIR

global globalFsDir;

if ~isempty(globalFsDir)
    fsurfSubDir=globalFsDir;
else
    if ispc,
        error('Hey mon, if you be using Windows you need to define the global variable "globalFsDir" and put the path to your FreeSurfer subject folder in it.');
    else
        fsurfSubDir=getenv('SUBJECTS_DIR');
        if isempty(fsurfSubDir)
            error('Could not read shell variable SUBJECTS_DIR. Try storing the path to your FreeSurfer subjects folder in the global variable globalFsDir.');
        end
    end
end


end

