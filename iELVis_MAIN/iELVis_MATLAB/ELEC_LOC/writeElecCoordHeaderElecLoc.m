function [fidXyz]=writeElecCoordHeader(fname)
%function [fidXyz]=writeElecCoordHeader(fname)
%
% Opens the file "fname" for writing and writes the header for electrode
% coordinates.
%
% Returns the file id, fidXyz, so that other functions can add coordinates
% to the file

fidXyz=fopen(fname,'w');
fprintf(fidXyz,'%s',datestr(now));
fprintf(fidXyz,'\n');
fprintf(fidXyz,'R A S\n');

% Future BIDS fields
%         fidXyz=fopen(fname,'w');
%     fprintf(fidXyz,'%s',datestr(now));
%     fprintf(fidXyz,'\tyangWang'); % method for brain shift correction
%     fprintf(fidXyz,'\tyangWang'); % version of freesurfer used
%     fprintf(fidXyz,'\n');
%     fprintf(fidXyz,'R A S\n');
end