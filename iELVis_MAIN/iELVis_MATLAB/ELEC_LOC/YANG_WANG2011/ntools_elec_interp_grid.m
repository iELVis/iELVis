function locpos = ntools_elec_interp_grid(ini_loc,ini_pos,row,col)

% use any 3 points in the grid that has a right angle to calculate the
% whole grid position. ini is the initial location of the 3 elecs, 3*3 
% matrix, in x y z cordinates. pos is the coresponding label for the 3 
% elecs in the grid, 3*1 vector. row and col are the numbers of the grid
% rows and columns. 

%Output:
% locpos = [pos,loc];
% pos: [row col] of the elecments in ascending order. 
% loc: xyz locations of the elements of the coresponding pos.


[mrow, mcol] = meshgrid(1:row,1:col);
pos_a = [mrow(ini_pos(1)),mcol(ini_pos(1))];
pos_b = [mrow(ini_pos(2)),mcol(ini_pos(2))];
pos_c = [mrow(ini_pos(3)),mcol(ini_pos(3))];

loc_a = ini_loc(1,:); % location for A
loc_b = ini_loc(2,:); % location for B
loc_c = ini_loc(3,:); % location for C

%% Calculate the 4 edge points, seperate the grid into 4 parts
t = find(pos_a-pos_b);

%% If A & B are in the same row and B & C in the same column
if t==2 
% Points A & B are in the same row, D is the first one of the row, and E is
% the last one of the row.
    lamda_b = pos_b(2)-1; % how far from B to D
    lamda_a = pos_a(2)-1; % how far from A to D
    l = -lamda_b/lamda_a;
    if lamda_a==0 % D is same with A
        pos_d = pos_a;
        loc_d = loc_a;
    else
        pos_d = [pos_a(1),1];
        loc_d = (loc_b+l*loc_a)/(1+l);
    end

    lamda_b = col-pos_b(2); % how far from B to E
    lamda_a = col-pos_a(2); % how far from A to E
    l = -lamda_b/lamda_a;
    if lamda_a==0 % E is same with A
        pos_e = pos_a;
        loc_e = loc_a;
    else
        pos_e = [pos_a(1),col];
        loc_e = (loc_b+l*loc_a)/(1+l);
    end

% Points B & C are in the same column, F is the first one of the column, G
% is the last one of the column
    lamda_b = pos_b(1)-1; % how far from B to F
    lamda_c = pos_c(1)-1; % how far from C to F
    l = -lamda_b/lamda_c;
    if lamda_c==0 % F is same with C
        pos_f = pos_c;
        loc_f = loc_c;
    else
        pos_f = [1,pos_c(2)];
        loc_f = (loc_b+l*loc_c)/(1+l);
    end

    lamda_b = row-pos_b(1); % how far from B to G
    lamda_c = row-pos_c(1); % how far from C to G
    l = -lamda_b/lamda_c;
    if lamda_c==0 % G is same with C
        pos_g = pos_c;
        loc_g = loc_c;
    else
        pos_g = [row,pos_c(2)];
        loc_g = (loc_b+l*loc_c)/(1+l);
    end

%% If A & B are in the same col and B & C in the same row
else 
% Points A & B are in the same col, F is the first one of the col, and G is
% the last one of the col.
    lamda_b = pos_b(1)-1; % how far from B to F
    lamda_a = pos_a(1)-1; % how far from A to F
    l = -lamda_b/lamda_a;
    if lamda_a==0 % F is same with A
        pos_f = pos_a;
        loc_f = loc_a;
    else
        pos_f = [1,pos_a(2)];
        loc_f = (loc_b+l*loc_a)/(1+l);
    end

    lamda_b = row-pos_b(1); % how far from B to G
    lamda_a = row-pos_a(1); % how far from A to G
    l = -lamda_b/lamda_a;
    if lamda_a==0 % G is same with A
        pos_g = pos_a;
        loc_g = loc_a;
    else
        pos_g = [row,pos_a(2)];
        loc_g = (loc_b+l*loc_a)/(1+l);
    end

% Points B & C are in the same row, D is the first one of the column, E
% is the last one of the column
    lamda_b = pos_b(2)-1; % how far from B to D
    lamda_c = pos_c(2)-1; % how far from C to D
    l = -lamda_b/lamda_c;
    if lamda_c==0 % D is same with C
        pos_d = pos_c;
        loc_d = loc_c;
    else
        pos_d = [pos_c(1),1];
        loc_d = (loc_b+l*loc_c)/(1+l);
    end

    lamda_b = col-pos_b(2); % how far from B to E
    lamda_c = col-pos_c(2); % how far from C to E
    l = -lamda_b/lamda_c;
    if lamda_c==0 % E is same with C
        pos_e = pos_c;
        loc_e = loc_c;
    else
        pos_e = [pos_c(1),col];
        loc_e = (loc_b+l*loc_c)/(1+l);
    end
end

%% Do the ntools_elec_locpos 4 times and get the 4 grids' locations
ini_loc1 = [loc_e;loc_b;loc_g];
ini_pos1 = [pos_e;pos_b;pos_g];
[loc1, pos1] = ntools_elec_locpos(ini_loc1,ini_pos1);

ini_loc2 = [loc_e;loc_b;loc_f];
ini_pos2 = [pos_e;pos_b;pos_f];
[loc2, pos2] = ntools_elec_locpos(ini_loc2,ini_pos2);

ini_loc3 = [loc_d;loc_b;loc_f];
ini_pos3 = [pos_d;pos_b;pos_f];
[loc3, pos3] = ntools_elec_locpos(ini_loc3,ini_pos3);

ini_loc4 = [loc_d;loc_b;loc_g];
ini_pos4 = [pos_d;pos_b;pos_g];
[loc4, pos4] = ntools_elec_locpos(ini_loc4,ini_pos4);

loc5 = [loc1;loc2;loc3;loc4];
[loc, m] = unique(loc5,'rows');
pos5 = [pos1;pos2;pos3;pos4];
for i=1:length(pos5)
    num(i) = find(mrow==pos5(i,1) & mcol==pos5(i,2));
end
pos = num(m(:))';

% pos: [row col] of the elecments in ascending order. 
% loc: xyz locations of the elements of the coresponding pos.
locpos = [pos,loc];
locpos = sortrows(locpos);
