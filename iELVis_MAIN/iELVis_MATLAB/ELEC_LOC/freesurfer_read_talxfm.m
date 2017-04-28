function [TalairachXFM] = freesurfer_read_talxfm(file)
% freesurfer_read_talxfm - Read FreeSurfer talairach.xfm file
% 
% [TalairachXFM] = freesurfer_read_talxfm(file)
% 
% This function will read an ascii file that contains the talairach
% transformation matrix for a FreeSurfer MRI volume.  These files are
% located in <subject>/mri/transforms/talairach.xfm
%
% The transformation matrix is based on an affine coregistration of the
% FreeSurfer MRI volume with the MNI template volume.  To use the
% TalairachXFM with a surface file, see freesurfer_surf2tal or just
% consider the following:
%
% Lets assume the surface RAS coordinates are in a vector,
% v = [R A S 1]', and the talairach.xfm matrix is X, then MNI Talairach 
% coordinates are:
%
% vMNI = X * v;
%

% $Revision: 1.1 $ $Date: 2004/11/17 21:03:27 $

% Licence:  GNU GPL, no implied or express warranties
% History:  11/2004 Darren.Weber_at_radiology.ucsf.edu
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fid = fopen(file,'r');

if isequal(fid,-1),
    S=sprintf('Could not open file: "%s"',file);
    error(S);
else
    
    %fprintf('...Reading FreeSurfer talairach.xfm file:\n%s\n',file);
    
    % read lines until we get the string 'Linear_Transform', which precedes
    % the data transformation matrix
    gotit = 0;
    string2match = 'Linear_Transform';
    for i=1:20, % read up to 20 lines, no more
        temp = fgetl(fid);
        if strmatch(string2match,temp),
            % we have the right line, so don't read any more
            gotit = 1;
            break;
        end
    end
    
    if gotit,
        % Read the transformation matrix (3x4).
        TalairachXFM = fscanf(fid,'%f',[4,3])';
        fclose(fid);
    else
        fclose(fid);
        error('failed to find ''Linear_Transform'' string in first 20 lines of xfm file.');
    end
    
end

return