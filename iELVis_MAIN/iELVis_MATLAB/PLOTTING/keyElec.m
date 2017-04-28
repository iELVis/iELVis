function yesno=keyElec(elec_name,elec_names)
%function yesno=keyElec(elec_name,elec_names)
%
% Indicates if an electrode is a grid corner of the first or last
% electrode in a strip.
% Right now it only works for 64 channel grids. 
%
% Sloppily written by David Groppe
% Mehtalab 2013

num_ids=find((elec_name>=48).*(elec_name<=57));
elec_num=str2num(elec_name(num_ids));

yesno=1;
if elec_num==1,
    yesno=1;
elseif strcmpi(elec_name,'G8') || strcmpi(elec_name,'G57')  || strcmpi(elec_name,'G64')
    yesno=1;
else
    n_elec=length(elec_names);
    for a=1:n_elec,
        num_ids2=find((elec_names{a}>=48).*(elec_names{a}<=57));
        if num_ids2(1)==num_ids(1)
            if strcmpi(elec_name(1:num_ids(1)-1),elec_names{a}(1:num_ids2(1)-1))
                if str2num(elec_name(num_ids))<str2num(elec_names{a}(num_ids2))
                    yesno=0;
                    break;
                end
            end
        end
    end
end
