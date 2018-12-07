function [fidXyz]=writeElecCoordHeader(fname)
    fidXyz=fopen(fname,'w');
    fprintf(fidXyz,'%s',datestr(now));
    fprintf(fidXyz,'\n');
    fprintf(fidXyz,'R A S\n');
%         fidXyz=fopen(fname,'w');
%     fprintf(fidXyz,'%s',datestr(now));
%     fprintf(fidXyz,'\tyangWang'); % method for brain shift correction
%     fprintf(fidXyz,'\tyangWang'); % version of freesurfer used
%     fprintf(fidXyz,'\n');
%     fprintf(fidXyz,'R A S\n');
end

% ?? TODO need to save FID when used in yangWang