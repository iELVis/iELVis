function [MNIv,TALv] = freesurfer_surf2tal(SURFv,TalairachXFM)
% freesurfer_surf2tal - convert freesurfer RAS to Talairach coordinates
%
% [MNIv,TALv] = freesurfer_surf2tal(SURFv,TalairachXFM)
%
% The input 'SURFv' are from freesurfer_read_surf; they must be Nx3
% vertex coordinates in freesurfer RAS coordinates. The second
% input is the matrix from the talairach.xfm file (see
% freesurfer_read_talxfm).
%
% This function converts the RAS coordinates into the MNI Talairach
% coordinates (MNIv).  It can also return the adjusted Talairach
% coordinates (TALv), based on the notes by Matthew Brett. This function
% calls mni2tal to apply the conversion from MNI template space to 'true'
% Talairach space.  See Matthew Brett's online discussion of this topic,
% http://www.mrc-cbu.cam.ac.uk/Imaging/Common/mnispace.shtml
%
% The surface RAS coordinates are arranged into a column vector, 
% SURFv = [R A S 1]', which is multiplied by the talairach.xfm
% matrix, TalairachXFM, to give the MNI Talairach coordinates:
%
% MNIv = TalairachXFM * SURFv;
% TALv = mni2tal(MNIv);
%

% $Revision: 1.1 $ $Date: 2004/11/17 21:03:27 $

% Licence:  GNU GPL, no express or implied warranties
% History:  11/2004, Darren.Weber_at_radiology.ucsf.edu
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%--------------------------------------------------------------------------
% transpose SURFv and add ones to the matrix

Nvertices = size(SURFv,1);

SURFv = SURFv';
SURFv(4,:) = ones(1, Nvertices);

%--------------------------------------------------------------------------
% Convert FreeSurfer RAS to Talairach coordinates

ver = '$Revision: 1.1 $';
fprintf('FREESURFER_SURF2TAL [v %s]\n',ver(11:15));

MNIv = TalairachXFM * SURFv;
MNIv = MNIv'; % return Nx3 matrix

if nargout > 1,
    TALv = mni2tal(MNIv);
end

return

