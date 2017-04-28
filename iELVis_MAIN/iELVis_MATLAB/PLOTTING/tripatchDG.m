function handle=tripatchDG(struct, nofigure, varargin)
%function handle=tripatchDG(struct, nofigure, varargin)
%
% A *very* slightly modified version of the freesurfer function tripatch.m

if nargin<2 | isempty(nofigure)
   figure
end
if nargin<3
   handle=trisurf(struct.tri, struct.vert(:, 1), struct.vert(:, 2), struct.vert(:, 3));
else
   if isnumeric(varargin{1})
      col=varargin{1};
      varargin(1)=[];
      if [1 3]==sort(size(col))
         col=repmat(col(:)', [size(struct.vert, 1) 1]);
      end
        handle=trisurf(struct.tri, struct.vert(:, 1), struct.vert(:, 2), struct.vert(:, 3), ...
         'FaceVertexCData', col, varargin{:});
      if length(col)==size(struct.vert, 1)
         set(handle, 'FaceColor', 'interp');
      end
   else
      handle=trisurf(struct.tri, struct.vert(:, 1), struct.vert(:, 2), struct.vert(:, 3), varargin{:});
   end
end
axis tight
axis equal
hold on
%David Groppe commented the following lines out because the cameratoolbar
%was making it impossible to make the surface interactive (besides rotating
%it)
% if version('-release')>=12 DG ??
%    cameratoolbar('setmode', 'orbit') <-what does that do?
% else 
rotate3d on
% end

% Make brain transparent
%set(handle,'facealpha',.25);
%set(handle,'edgealpha',.25);