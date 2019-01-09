function fsurf2BIDS(fnameStems,sourceDir,destDir,optional)
% function fsurf2BIDS(fnameStems,sourceDir,destDir,optional)
%
% Copies both hemisphere versions of a FreeSurfer file (e.g., lh.pial and 
% rh.pial) from sourceDir into destDir
%
% Required Inputs:
%  fnameStems - The name of FreeSurfer minus the lh/rh prefix (e.g., '.pial')
%  sourceDir  - The directory from which to get the file
%  destDir    - The destination directory
%
% Optional Input:
%  optional  - If non-zero, no errors will be thrown if file copying failes
%              (e.g., the targetted file doesn't exist). {default: 0}

if nargin<4,
   optional=0; 
end

[SUCCESS,MESSAGE,MESSAGEID] = mkdir(destDir);
for a=1:length(fnameStems),
    for b=1:2,
        if b==1,
            hem='lh';
        else
            hem='rh';
        end
        myFname=[hem '.' fnameStems{a}];
        if universalYes(optional),
            % won't throw an error
            try
                copyfile(fullfile(sourceDir,myFname),fullfile(destDir,myFname));
            catch
                
            end
        else
            % Throws an error if unsuccesful
            copyfile(fullfile(sourceDir,myFname),fullfile(destDir,myFname));
        end
    end
end

end