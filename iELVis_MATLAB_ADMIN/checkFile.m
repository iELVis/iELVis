function checkFile(fname)
%function checkFile(fname)

if ~exist(fname,'file')
   error('File %s not found.',fname); 
end