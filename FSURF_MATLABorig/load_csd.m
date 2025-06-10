function [dat thresh fwhm searchspace anattype] = load_csd(csdfile)
% [dat thresh fwhm searchspace anattype] = load_csd(csdfile)
% 
% reads in Cluster Simulation Data (CSD) as produced by mri_glmfit.
% dat has 5 columes:
%  1. Row number (0-based)
%  2. nClusters 
%  3. MaxClustSize
%  4. MaxSig    
%  5. MaxStat
% 
% thresh is -log10(p)
% fwhm in mm^D
% searchspace in mm^D
% anattype is volume or surface
%
% $Id: load_csd.m,v 1.3 2009/01/17 01:19:51 greve Exp $

%
% load_csd.m
%
% Original Author: Doug Greve
% CVS Revision Info:
%    $Author: greve $
%    $Date: 2009/01/17 01:19:51 $
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

if(nargin ~= 1)
  nargin
  fprintf('dat = load_csd(csdfile)\n');
  return;
end


%'synth-glm-surf/n0001/csd/mc-z.pos.j001-osgm.csd'
fid = fopen(csdfile,'r');
if(fid == -1)
  fprintf('ERROR: opening %s\n',csdfile);
  return;
end

tline = fgetl(fid);
if(tline == -1)
  fprintf('ERROR: %s is not correctly formatted, no first line\n', ...
	  csdfile);
  fclose(fid);
  return;
end

%----------- Loop through all the lines ----------------------%
nthrow = 1;
while(1)

  % scroll through any blank lines or comments %
  while(1)
    tline = fgetl(fid);
    if(~isempty(tline))
      if(tline(1) ~= '#') break; end
      key = sscanf(tline,'%*s %s',1);
      if(strcmp(key,'thresh'))
	thresh = sscanf(tline,'%*s %*s %f');
      end
      if(strcmp(key,'nullfwhm'))
	fwhm = sscanf(tline,'%*s %*s %f');
      end
      if(strcmp(key,'searchspace'))
	searchspace = sscanf(tline,'%*s %*s %f');
      end
      if(strcmp(key,'anattype'))
	anattype = sscanf(tline,'%*s %*s %s');
      end
    end
  end
  if(tline(1) == -1) break; end
  dat(nthrow,:) = sscanf(tline,'%f',5);

  nthrow = nthrow + 1;
end % while (1)

fclose(fid);

return;








