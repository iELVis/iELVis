function [fidXyz]=writeElecCoordHeader(fname,brainShiftMethod)
%function [fidXyz]=writeElecCoordHeader(fname)
%
% Inputs:
%  fname - the filename to be opened for writing
%  brainShiftMethod - a string that specifies the method used for
%             correcting for brain shift and the version of the code
%             (represented by the git hash string)
%
% Opens the file "fname" for writing and writes the header for electrode
% coordinates.
%
% Returns the file id, fidXyz, so that other functions can add coordinates
% to the file

% Orig info
% fidXyz=fopen(fname,'w');
% fprintf(fidXyz,'%s',datestr(now));
% fprintf(fidXyz,'\n');
% fprintf(fidXyz,'R A S\n');

fidXyz=fopen(fname,'w');
fprintf(fidXyz,'%s',datestr(now));
fprintf(fidXyz,'\t%s',brainShiftMethod); % method for brain shift correction
fprintf(fidXyz,'\n');
fprintf(fidXyz,'R A S\n');
% TODO add version of freesurfer used
end