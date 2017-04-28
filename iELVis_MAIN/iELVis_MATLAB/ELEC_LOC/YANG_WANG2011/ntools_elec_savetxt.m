function filename = ntools_elec_savetxt(filename,varargin)
% save the electrodes locations into a text file
%
% usage: ntools_elec_savetxt(output_filename, data1,data2,....)
% the class of data could be struct, cell or double
%
% the first column of data is the electrode name, the other 3
% columns are the x y z coordinates.
%
% if the data is a double matirx, then the input data name will be the 
% electrode name and be output in the first column of the final text file.

fid = fopen(filename,'w');
% fclose(fid);
% fid = fopen(filename,'a');
for k = 2:nargin
    elec = varargin{k-1};
    if ~isempty(elec)
        c = class(elec);
        switch c
            case 'struct'
                elecname = fieldnames(elec);
                for i = 1:length(elecname)
                    for j = 1:size(elec.(char(elecname{i})),1)
                        name_num = sprintf('%s%.2d',elecname{i},j);
                        fprintf(fid,'%s %f %f %f \n',upper(name_num),...
                            elec.(char(elecname{i}))(j,1),elec.(char(elecname{i}))(j,2),elec.(char(elecname{i}))(j,3));
                    end
                end
            case 'cell'
                for i = 1:size(elec,1)
                    %                 strip_name = regexp(elec(i,1),'[A-Za-z]*[^\d*]','match');
                    %                 strip_num = regexp(elec(i,1),'[^A-Za-z]*[\d*]','match');
                    %                 name_num = sprintf('%s%.2d',upper(char(strip_name{1})),str2double(strip_num{1}));
                    %                 fprintf(fid,'%s %f %f %f \n',name_num,elec{i,2},elec{i,3},elec{i,4});
                    fprintf(fid,'%s %f %f %f %s \n',elec{i,1},elec{i,2},elec{i,3},elec{i,4},elec{i,5});
                end
            case 'double'
                elecname = upper(inputname(k));
                %             name_num = sprintf('%s%.2d',elecname,i);
                for i = 1:size(elec,1)
                    name_num = sprintf('%s%.2d',elecname,i);
                    fprintf(fid,'%s %f %f %f\n',name_num,elec(i,1),elec(i,2),elec(i,3));
                end
        end
    end

    clear elec;

end
fclose(fid);
% [name x y z] = textread(filename,'%s %f %f %f');
% electrodes = [x y z];

return;