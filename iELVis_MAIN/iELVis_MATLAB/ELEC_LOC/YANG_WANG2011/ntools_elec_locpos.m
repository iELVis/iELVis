function [loc, pos] = ntools_elec_locpos(ini_loc,ini_pos)

% ------------------------------------
% This program is used to calculate the position of 
% every matrix element in a 3D or 2D plane or straight line
% 
%           A     col     B
%             *   *   *   *                             col
%  row        *   *   *   *        or       A * * * * * * * * B
%             *   *   *   *
%             *   *   *   *
%            D             C
% 
% Input:
% ini_loc: the initial points xyz locations. for a plane, it must be
%               [loc_A;loc_B;loc_C], a 3*3 or 3*2 matrix; for a straight line, it must be
%               [loc_A;loc_B], a 2*3 or 2*2 matrix.
% ini_pos: the row and column numbers of the A;B;C in the grid, 3*2 matrix. 
%               for the straight line, it's the index number of A & B, 2*1 matrix;             
% 
% Output:
% pos: [row col] of the elecments in ascending order. 
% loc: xyz locations of the elements of the coresponding pos.
% ------------------------------------

% Examples
% ini_loc = [-47.6667  -15.3333   47.9967; -29.1000  -58.9000   39.4300; -38.4200  -59.5800   21.4300];
% ini_loc = [-15.3333   47.9967; -58.9000   39.4300; -59.5800   21.4300];
% ini_pos = [2,1;2,6; 4,6];
%
% ini_loc = [-40.564110 49.700199 17.898691;-45.409290 -12.858190 52.881947;-11.667823 -48.278912 -10.167865];
% ini_pos = [1,1;1,8;8,8];
%
% ini_loc = [170  144  116; 150  123   114];
% ini_loc = [-15.3333   47.9967; -58.9000   39.4300];
% ini_pos = [8;1];


form= size(ini_loc);
D3 = [3,3]; D2 = [3,2]; D1 = [2,3]; D0 = [2,2];

if isequal(form,D3) || isequal(form,D2) % 3D or 2D grid interp
    
    if ~isequal(size(ini_pos),D2)
        error('the initial position matrix should be a 3*2 matrix');
    end
    
    loc_a = ini_loc(1,:); % A
    loc_b = ini_loc(2,:); % B
    loc_c = ini_loc(3,:); % C
    pos_a = ini_pos(1,:);
    pos_b = ini_pos(2,:);
    pos_c = ini_pos(3,:);


    %% Check if it construct a matix
    if isequal(pos_a,pos_b) || isequal(pos_c,pos_b)
         loc = []; pos = []; return;
    end

    % Calculate the rows and cols
    col = abs(pos_a(2)-pos_b(2))+1; % a & b are in the same row
    row = abs(pos_c(1)-pos_b(1))+1; % b & c are in the same column

    %% Determine the 4th point D which has the same row with C, and same
    %% column with A
    BA = loc_a - loc_b;
    BC = loc_c - loc_b;
    BD = BA + BC;
    loc_d = BD + loc_b;
%     pos_d = [pos_c(1),pos_a(2)];

    %% Use interp2 for x,y,z seperately
    r_max = max([pos_a(1);pos_b(1);pos_c(1)]);
    r_min = min([pos_a(1);pos_b(1);pos_c(1)]);
    c_max = max([pos_a(2);pos_b(2);pos_c(2)]);
    c_min = min([pos_a(2);pos_b(2);pos_c(2)]);

    [x, y] = meshgrid(c_min:c_max,r_min:r_max); 
    [X, Y] = meshgrid([pos_a(2),pos_c(2)],[pos_a(1),pos_c(1)]);
    
    if form(2)==3
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
        
    elseif form(2)==2
        Zx = [loc_a(1),loc_b(1);loc_d(1),loc_c(1)];
        Zy = [loc_a(2),loc_b(2);loc_d(2),loc_c(2)];
        
        loc_x = interp2(X,Y,Zx,x,y);
        loc_y = interp2(X,Y,Zy,x,y);
        
        loc = zeros(row*col,2);
        pos = zeros(row*col,2);
        loc_temp = zeros([size(loc_x),2]);
        loc_temp(:,:,1) = loc_x;
        loc_temp(:,:,2) = loc_y;
    end
    
    k = 1;
    for i=1:row
        for j=1:col
            loc(k,:) = loc_temp(i,j,:);
            pos(k,:) = [y(i,j),x(i,j)];
            k = k+1;
        end
    end
    
elseif isequal(form,D1) || isequal(form,D0) % 3D or 2D straight line interp
   
    if ~isequal(size(ini_pos),[2 1])
        error('the initial position matrix should be a 2*1 matrix');
    end

    col = abs(ini_pos(1)-ini_pos(2))+1;
    
    %% Check if it construct a matix
    if col==1
        loc = []; pos = []; return;
    end
    
    %% use interp1 for straight line interp
    c_max = max(ini_pos);
    c_min = min(ini_pos);
    pos = (c_min:c_max)';
    if form(2)==3
        loc_x = interp1(ini_pos,ini_loc(:,1),pos);
        loc_y = interp1(ini_pos,ini_loc(:,2),pos);
        loc_z = interp1(ini_pos,ini_loc(:,3),pos);
        loc = [loc_x,loc_y,loc_z];
    elseif form(2)==2
        loc_x = interp1(ini_pos,ini_loc(:,1),pos);
        loc_y = interp1(ini_pos,ini_loc(:,2),pos);
        loc = [loc_x,loc_y];
    end
    
else
    error('initial location must be one of a 3*3, 3*2, 2*3, 2*2 matrix');
end
