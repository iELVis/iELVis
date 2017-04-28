function locpos = ntools_elec_locpos4(ini_loc,ini_pos,row,col)

% Interpolate the grid electrodes locations with 4 initial elecs. The
% initial elecs have to be the corner ones in either clockwise or
% counter-clockwise direction
%
% The difference from ntools_elec_locpos:
% 1. only takes 4*3 initial elecs
% 2. only takes corner elecs
% 3. ini_pos are elec index, instead of [row, col] in ntools_elec_locpos
% 4. output is the same with ntools_elec_interp_grid
%
% created by Hugh Wang on Oct.8th, 2013

% convert from index to row,col
[mrow, mcol] = meshgrid(1:row,1:col);
pos_a = [mrow(ini_pos(1)),mcol(ini_pos(1))];
pos_b = [mrow(ini_pos(2)),mcol(ini_pos(2))];
pos_c = [mrow(ini_pos(3)),mcol(ini_pos(3))];
pos_d = [mrow(ini_pos(4)),mcol(ini_pos(4))];

loc_a = ini_loc(1,:);
loc_b = ini_loc(2,:);
loc_c = ini_loc(3,:);
loc_d = ini_loc(4,:);

%% Use interp2 for x,y,z seperately
r_max = max([pos_a(1);pos_b(1);pos_c(1);pos_d(1)]);
r_min = min([pos_a(1);pos_b(1);pos_c(1);pos_d(1)]);
c_max = max([pos_a(2);pos_b(2);pos_c(2);pos_d(2)]);
c_min = min([pos_a(2);pos_b(2);pos_c(2);pos_d(2)]);

[x, y] = meshgrid(c_min:c_max,r_min:r_max);
[X, Y] = meshgrid([pos_a(2),pos_c(2)],[pos_a(1),pos_c(1)]);

Zx = [loc_a(1),loc_b(1);loc_d(1),loc_c(1)];
Zy = [loc_a(2),loc_b(2);loc_d(2),loc_c(2)];
Zz = [loc_a(3),loc_b(3);loc_d(3),loc_c(3)];

loc_x = interp2(X,Y,Zx,x,y);
loc_y = interp2(X,Y,Zy,x,y);
loc_z = interp2(X,Y,Zz,x,y);

%% Save the data to loc & pos
loc = zeros(row*col,3);
pos = zeros(row*col,2);
loc_temp = zeros([size(loc_x),3]);
loc_temp(:,:,1) = loc_x;
loc_temp(:,:,2) = loc_y;
loc_temp(:,:,3) = loc_z;

k = 1;
for i=1:row
    for j=1:col
        loc(k,:) = loc_temp(i,j,:);
        pos(k,:) = [y(i,j),x(i,j)];
        k = k+1;
    end
end

[loc1, m] = unique(loc,'rows');
for i=1:length(pos)
    num(i) = find(mrow==pos(i,1) & mcol==pos(i,2));
end
pos1 = num(m(:))';

locpos = [pos1,loc1];
locpos = sortrows(locpos);

end

