function elec_names=formatElecNames(elec_names)
%function elec_names=formatElecNames(elec_names)
%
% Removes underscores AND spaces from electrode names and replaces 'Grid' with 'G'
% Input:
%   elec_names - cell array of electrode names
%
% Output:
%   elec_names - reformatted cell array of electrode names
%
% D. Groppe


n_chan=length(elec_names);
for a=1:n_chan,
    % Remove underscores
    n_char=length(elec_names{a});
    temp_str=repmat('a',1,n_char);
    ct=0;
    for b=1:n_char,
        if elec_names{a}(b)~='_' && elec_names{a}(b)~=' '
            ct=ct+1;
            temp_str(ct)=elec_names{a}(b);
        end
    end
    elec_names{a}=temp_str(1:ct);
    if length(elec_names{a})>4
        if strcmpi(elec_names{a}(1:4),'Grid')
            elec_names{a}=['G' elec_names{a}(5:end)];
        end
    end
end