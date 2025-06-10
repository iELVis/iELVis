

%
% talairaching_table.m
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


dirname='/space/thigpen/1/users/efenster/PPG_4_28_05';
tabledir='/tmp';
% Creates the tables in tabledir %
[D,mu,sigma]=talairaching_stats(dirname, tabledir); % stores TalairachingMean.adf and TalairachingCovariance.adf in the outdir (arg2)
mean=strcat(tabledir, '/TalairachingMean_tmp.adf');
sigma=strcat(tabledir, '/TalairachingCovariance_tmp.adf');
mean2=strcat(tabledir, '/TalairachingMean.adf');
sigma2=strcat(tabledir, '/TalairachingCovariance.adf');
fp_mu=fopen('/tmp/header_mu.txt', 'w');
fprintf(fp,'# FSAFD TalairachingCheck 1\n# date 20050603\n# $Id: talairaching_table.m,v 1.3 2007/01/10 22:55:10 nicks Exp $\n# efenster PPG_4_28_05 \n# nsubjects %d\n# ncols 9\n# nrows 1\n# Average talairach.xfm coefficients (translation excluded)\n', length(D));
fclose(fp_mu);
fp_sigma=fopen('/tmp/header_sigma.txt', 'w');
fprintf(fp,'# FSAFD TalairachingCheck 1\n# date 20050603\n# $Id: talairaching_table.m,v 1.3 2007/01/10 22:55:10 nicks Exp $\n# efenster PPG_4_28_05 \n# nsubjects %d\n# ncols 9\n# nrows 9\n# covariance matrix\n', length(D));
fclose(fp_mu);
unix(sprintf('cat /tmp/header_mu.txt %s > %s', mean, mean2));
unix(sprintf('cat /tmp/header_sigma.txt %s > %s', sigma, sigma2));
unix(sprintf('rm %s | rm %s', mean, sigma));
unix(sprintf('rm %s | rm %s', '/tmp/header_mu.txt', '/tmp/header_sigma.txt'));

[P, I, nf]=talairaching_dir_afd(dirname, 0.01, tabledir); % stores TalairachingProbas.adf in the outdir (arg3)
%load('/space/okapi/3/data/laurence/ADF/talairaching/transfo_param_probas.mat') %loads y computed for the Buckner data set
%P=y;
%load('/space/okapi/3/data/laurence/ADF/talairaching/transfo_param_names.mat') %loads nf (Buckner data set)
files=dir(dirname);
table_proba=strcat(tabledir, '/TalairachingProbas.adf');
fp=fopen(table_proba, 'w');
fprintf(fp,'# FSAFD TalairachingCheck 1\n# date 20050603\n# $Id: talairaching_table.m,v 1.3 2007/01/10 22:55:10 nicks Exp $\n# efenster PPG_4_28_05 \n# ncols 1\n# nrows %d\n# label_col 1 xfmProbability\n', length(P));
for i=1:length(P)
    %u=I(i);
    %fprintf(fp, '# label_row %d %s\n', i, char(files(u).name));
    fprintf(fp, '# label_row %d %s\n', i, char(nf(i).name));
end
for j=1:length(P)
    fprintf(fp, '%f\n', P(j));
end
fclose(fp)
