function hdr = ntools_elec_load_nifti(niftifile,hdronly)
% hdr = load_nifti(niftifile,hdronly)
%
% Loads nifti header and volume. The volume is stored
% in hdr.vol. Columns and rows are not swapped.
%
% Handles compressed nifti (nii.gz) by issuing a unix command to
% uncompress the file to a temporary file, which is then deleted.
%
% Dimensions are in mm and msec
% hdr.pixdim(1) = physical size of first dim (eg, 3.125 mm or 2000 ms)
% hdr.pixdim(2) = ...
% 
% The sform and qform matrices are stored in hdr.sform and hdr.qform.
%
% hdr.vox2ras is the vox2ras matrix based on sform (if valid), then
% qform.
%
% See also: load_nifti_hdr.m
%


%
% load_nifti.m
%
% Original Author: Doug Greve
% CVS Revision Info:
%    $Author: nicks $
%    $Date: 2007/12/10 15:52:36 $
%    $Revision: 1.13.2.1 $
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

hdr = [];

if(nargin < 1 || nargin > 2)
  fprintf('hdr = load_nifti(niftifile,<hdronly>)\n');
  return;
end

if(~exist('hdronly','var')) hdronly = []; end
if(isempty(hdronly)) hdronly = 0; end

% unzip if it is compressed 
ext = niftifile((strlen(niftifile)-2):strlen(niftifile));
if(strcmpi(ext,'.gz'))
  % Need to create unique file name (harder than it looks)
  rand('state', sum(100*clock));
  gzipped =  round(rand(1)*10000000 + ...
		   sum(int16(niftifile))) + round(cputime);
  ind = findstr(niftifile, '.');
  if strncmp(computer,'PCWIN',5) % PM edited 20160322
      new_niftifile = sprintf('tmp%d.nii', gzipped);
  else
      new_niftifile = sprintf('/tmp/tmp%d.nii', gzipped);
  end
  %fprintf('Uncompressing %s to %s\n',niftifile,new_niftifile);
  if strncmp(computer,'PCWIN',5) % PM edited 20160322
      if ~exist('C:\Program Files\7-Zip\7z.exe','file')
          error('Win:no7zip','Could not find 7-Zip. Make sure that it is installed in\n''C:\\Program Files\\7-Zip\\'', or else edit this function\nto direct it to the correct location on your machine.');
      end
      % get folder and file name
      [foldername,fnameshort]=fileparts(niftifile);
      % execute 7-Zip extract command
      eval(['dos ''"C:\Program Files\7-Zip\7z" e "' niftifile '" -o"' foldername '"''']);
      % rename unzipped file
      movefile(fullfile(foldername,fnameshort),fullfile(foldername,new_niftifile));
  elseif(strcmp(computer,'MAC') || strcmp(computer,'MACI64') || strcmp(computer,'MACI'))
    unix(sprintf('gunzip -c %s > %s', niftifile, new_niftifile));
  else
    unix(sprintf('zcat %s > %s', niftifile, new_niftifile)) ;
  end
  niftifile = new_niftifile ;
else
  gzipped = -1 ;
end

if strncmp(computer,'PCWIN',5) % PM edited 20160322
    hdr=load_nifti_hdr(fullfile(foldername,niftifile));
else
    hdr = load_nifti_hdr(niftifile);
end
if(isempty(hdr)) 
  if(gzipped >=0)
      if strncmp(computer,'PCWIN',5) % PM edited 20160322
          delete(fullfile(foldername,new_niftifile));
      else
          unix(sprintf('rm %s', niftifile));
      end
  end
  return; 
end

% If only header is desired, return now
if(hdronly) 
  if(gzipped >=0)
      if strncmp(computer,'PCWIN',5) % PM edited 20160322
          delete(fullfile(foldername,new_niftifile));
      else
          unix(sprintf('rm %s', niftifile));
      end
  end
  return; 
end

% Open to read the pixel data
fp = fopen(niftifile,'r',hdr.endian);

% Get past the header
fseek(fp,round(hdr.vox_offset),'bof');

switch(hdr.datatype)
 % Note: 'char' seems to work upto matlab 7.1, but 'uchar' needed
 % for 7.2 and higher. 
 case   2, [hdr.vol, nitemsread] = fread(fp,inf,'uchar');
 case   4, [hdr.vol, nitemsread] = fread(fp,inf,'short');
 case   8, [hdr.vol, nitemsread] = fread(fp,inf,'int');
 case  16, [hdr.vol, nitemsread] = fread(fp,inf,'float');
 case  64, [hdr.vol, nitemsread] = fread(fp,inf,'double');
 case 512, [hdr.vol, nitemsread] = fread(fp,inf,'ushort');
 case 768, [hdr.vol, nitemsread] = fread(fp,inf,'uint');
 otherwise,
  fprintf('ERROR: data type %d not supported',hdr.datatype);
  hdr = [];
  return;
end

fclose(fp);
if(gzipped >=0)
    if strncmp(computer,'PCWIN',5) % PM edited 20160322
        delete(fullfile(foldername,new_niftifile));
    else
        unix(sprintf('rm %s', niftifile));
    end
end

% Get total number of voxels
dim = hdr.dim(2:end);
ind0 = dim==0;
dim(ind0) = 1;
nvoxels = prod(dim);

% Check that that many voxels were read in
if(nitemsread ~= nvoxels) 
  fprintf('ERROR: %s, read in %d voxels, expected %d\n',...
	  niftifile,nitemsread,nvoxels);
  hdr = [];
  return;
end

hdr.vol = reshape(hdr.vol, dim');
if(hdr.scl_slope ~= 0)
  fprintf('Rescaling NIFTI: slope = %g, intercept = %g\n',...
	  hdr.scl_slope,hdr.scl_inter);
  hdr.vol = hdr.vol * hdr.scl_slope  + hdr.scl_inter;
end

return;