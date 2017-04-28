function tpt=findTpt(tme,tmes)
%function tpt=findTpt(tme,tmes)
%
% Inputs:
%   tme  - [scalar] a desired time
%   tmes - a vector of times
%
% Output:
%   tpt  - the element of tmes that is closest to tme in terms of absolute
%          difference (L1 norm)
%
% Note: tme and tmes should be in the same units (e.g., milliseconds)

[dummy tpt]=min(abs(tme-tmes));
