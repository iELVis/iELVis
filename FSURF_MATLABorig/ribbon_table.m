

%
% ribbon_table.m
%
% Original Author: Laurence Wastiaux
% CVS Revision Info:
%    $Author: nicks $
%    $Date: 2007/01/10 22:55:10 $
%    $Revision: 1.2 $
%
% Copyright (C) 2002-2007,
% The General Hospital Corporation (Boston, MA). 
% All rights reserved.
%
% Distribution, usage and copying of this software is covered under the
% terms found in the License Agreement file named 'COPYING' found in the
% FreeSurfer source code root directory, and duplicated here:
% https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferOpenSourceLicense
%
% General inquiries: freesurfer@nmr.mgh.harvard.edu
% Bug reports: analysis-bugs@nmr.mgh.harvard.edu
%



[D, Isubj]=ribbon_dir_afd('/space/neo/2/recon/buckner', 0.01);
%load('/space/okapi/3/data/laurence/ADF/ribbon/Dice_ribbon.mat') %loads D
%load('/space/okapi/3/data/laurence/ADF/ribbon/Subj_name_ribbon.mat') %loads nf
files=dir('/space/neo/2/recon/buckner');
fp=fopen('ribbon_table.txt', 'w');
fprintf(fp,'# FSAFD RibbonPlacementCheck 1\n# date 20050516\n# $Id: ribbon_table.m,v 1.2 2007/01/10 22:55:10 nicks Exp $\n# Buckner data set\n# ncols 1\n# nrows %d\n# label_col 1 DiceCoefficient\n', length(D));
for i=1:length(D)
    %nb=Isubj(i);
    %fprintf(fp, '# label_row %d %s\n', i, char(files(nb).name));
    fprintf(fp, '# label_row %d %s\n', i, char(nf(i).name));
end
for j=1:length(D)
    fprintf(fp, '%f\n', D(j));
end
fclose(fp)
