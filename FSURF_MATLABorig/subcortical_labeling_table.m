

%
% subcortical_labeling_table.m
%
% Original Author: Laurence Wastiaux
% CVS Revision Info:
%    $Author: nicks $
%    $Date: 2007/01/10 22:55:10 $
%    $Revision: 1.3 $
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



[Dvolume,Ivol]=subcortical_labeling_dir_afd('/space/neo/2/recon/buckner', 0.01)
% load('/space/okapi/3/data/laurence/ADF/subcortical_labeling/Percent_Vol0517.mat') %loads Dvolume
% load('/space/okapi/3/data/laurence/ADF/subcortical_labeling/Subj_nb_Vol0517.mat') %loads Ivol
files=dir('/space/neo/2/recon/buckner');
fp=fopen('~/subcortlab_table.txt', 'w');
sz=size(Dvolume);
fprintf(fp,'# FSAFD SubcorticalLabelingCheck 1\n# date 20050516\n# $Id: subcortical_labeling_table.m,v 1.3 2007/01/10 22:55:10 nicks Exp $\n# Buckner data set\n# ncols 20\n# nrows %d\n', sz(1));
ROI_Labels=[4 43 17 53 10 49 11 50 12 51 13 52 18 54 26 58 14 15 5 44];
labmapfile=('/space/lyon/1/fsdev/freesurfer_dev');
labmap=strcat(labmapfile,'/','tkmeditColorsCMA');
[label name val1 val2 val3 val4]=textread(labmap,'%d %s %d %d %d %d',89);
nn=1;
for u=1:length(ROI_Labels)
    fprintf(fp, '# label_col %d %s\n', u, char(name(ROI_Labels(u)+1)));
end
for i=1:sz(1)
     if(Dvolume(i,1)>0)
        nb=Ivol(i);
        fprintf(fp, '# label_row %d %s\n',nn, char(files(nb).name));
         nn=nn+1;
     end
end
for j=1:sz(1)
    for k=1:20
        fprintf(fp,'%f  ', Dvolume(j,k));
    end
    fprintf(fp, '\n');
end
fclose(fp)
