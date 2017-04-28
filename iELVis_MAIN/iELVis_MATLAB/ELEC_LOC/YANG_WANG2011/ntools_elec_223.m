function [loc,pos] = ntools_elec_223(ini_loc,ini_pos,row,col,sph)

% Calculate the elecs from 2 points to 3 points. ini_loc is a 2*3 vector,
% ini_pos is a 2*1 vector with position numbers. basically it solves 2
% equations: 1) dot product of the two vectors are 0 since the right angle,
% and 2) the norms of the vectors proportional which the ratio can be
% learnt from the ini_pos
%
% if the 2 initial points are not in the same row and not in the same column,  
% Output B (the right angle point) will be in the same row with A and same column with C. 
%
% if the 2 initial points are in the same row or in the same column, the
% output B has the position pos_c+1, C is the right angle point

% test initial
% clear all; close all;
% row=4; col=8;
% sph = 'lh';
% 
% % ini_pos = [11;22];
% % ini_loc = [-40.24 -32.76 44.57; -33.76 -59.24 30.43]; 
% 
% ini_pos = [1;2];
% ini_loc = [-40.24 -32.76 44.57; -37 -60.834 39.788];

%% Calculate the 3rd point
[mrow, mcol] = meshgrid(1:row,1:col);
pos_a = [mrow(ini_pos(1)),mcol(ini_pos(1))];
pos_c = [mrow(ini_pos(2)),mcol(ini_pos(2))];
dif = pos_a-pos_c;
if dif(1)~=0 && dif(2)~=0
    pos_b = [pos_a(1),pos_c(2)];
    num = find(mrow==pos_b(1)& mcol==pos_b(2));
    pos = [ini_pos(1);num;ini_pos(2)];
    label = 1;
elseif dif(1)==0 
    if pos_c(1)==row
        pos_b = [pos_c(1)-1,pos_c(2)];
    else
        pos_b = [pos_c(1)+1,pos_c(2)];
    end
    num = find(mrow==pos_b(1)& mcol==pos_b(2));
    pos = [ini_pos(1);ini_pos(2);num];
    label = 0;
elseif dif(2)==0
    if pos_c(2)==col
        pos_b = [pos_c(1),pos_c(2)-1];
    else
        pos_b = [pos_c(1),pos_c(2)+1];
    end 
    num = find(mrow==pos_b(1)& mcol==pos_b(2));
    pos = [ini_pos(1);ini_pos(2);num];
    label = 0;
end
%% Solve the equations
loc_a = ini_loc(1,:); % location for A
loc_c = ini_loc(2,:); % location for C

x = (loc_a(1)+loc_c(1))/2; % B (x)
y1 = num2str(loc_a(2)); % A (y)
z1 = num2str(loc_a(3)); % A (z)
y2 = num2str(loc_c(2)); % C (y)
z2 = num2str(loc_c(3)); % C (z)
if label==1
    bc2 = num2str(sum(pos_b-pos_c)^2); % bc
    ba2 = num2str(sum(pos_b-pos_a)^2); % ba
    a = sprintf('(%s-y)*(%s-y)+(%s-z)*(%s-z)',y1,y2,z1,z2); 
    b = sprintf('%s*(%s-y)^2+%s*(%s-z)^2-%s*(%s-y)^2-%s*(%s-z)^2',bc2,y1,bc2,z1,ba2,y2,ba2,z2);
elseif label==0
    cb2 = num2str(sum(pos_c-pos_b)^2); % cb
    ca2 = num2str(sum(pos_c-pos_a)^2); % ca
    a = sprintf('(y-(%s))*((%s)-(%s))+(z-(%s))*((%s)-(%s))',y2,y1,y2,z2,z1,z2);
    b = sprintf('%s*(y-(%s))^2+%s*(z-(%s))^2-%s*((%s)-(%s))^2-%s*((%s)-(%s))^2',...
        ca2,y2,ca2,z2,cb2,y1,y2,cb2,z1,z2);
end
s = solve(sym(a),sym(b),'y','z');
y = double(s.y);
z = double(s.z);

%% Plot the points
cmd1 = sprintf('Grid # %s',num2str(ini_pos(1)));
cmd2 = sprintf('Grid # %s',num2str(num));
cmd3 = sprintf('Grid # %s',num2str(ini_pos(2)));
plot([ini_loc(1,2);y(1);ini_loc(2,2);ini_loc(1,2)],[ini_loc(1,3);z(1);ini_loc(2,3);ini_loc(1,3)],...
    '--b','LineWidth',2);axis equal;
text(ini_loc(1,2)+1,ini_loc(1,3)+1,cmd1,'FontSize',12);
text(y(1)+1,z(1)+1,cmd2,'FontSize',12);
text(ini_loc(2,2)+1,ini_loc(2,3)+1,cmd3,'FontSize',12);

if strcmp(sph,'lh')
    set(gca,'XDir','reverse');
end
a = menu('Does the point you want locates like this?','Yes','No');
if a==1
    if label==1
        loc = [loc_a;[x,y(1),z(1)];loc_c];
    elseif label==0
        loc = [loc_a;loc_c;[x,y(1),z(1)]];
    end
else
    plot([ini_loc(1,2);y(2);ini_loc(2,2);ini_loc(1,2)],[ini_loc(1,3);z(2);ini_loc(2,3);ini_loc(1,3)],...
        '--b','LineWidth',2);axis equal;
    text(ini_loc(1,2)+1,ini_loc(1,3)+1,cmd1,'FontSize',12);
    text(y(2)+1,z(2)+1,cmd2,'FontSize',12);
    text(ini_loc(2,2)+1,ini_loc(2,3)+1,cmd3,'FontSize',12);
    if strcmp(sph,'lh')
        set(gca,'XDir','reverse');
    end
    b = menu('Does the point you want locates like this?','Yes','No');
    if b==1
        if label==1
            loc = [loc_a;[x,y(2),z(2)];loc_c];
        elseif label==0
            loc = [loc_a;loc_c;[x,y(2),z(2)]];
        end
    else 
        disp('Sorry, none satisfies your need. Please re-do the process.');
        loc = []; pos = [];
    end
end
close all;
return    
