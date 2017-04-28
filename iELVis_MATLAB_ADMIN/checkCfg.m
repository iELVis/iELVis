% check_cfg()-Throws an error if a function is being called with a
%             cfg-based input argument that is not defined by the function.
%
% Usage:
%  >>check_cfg(cfg,mfile);
%
% Required Inputs:
%  cfg   - [struct variable] The cfg variable with fields corresponding to 
%          input parameters
%  mfile - [string] The mfile that will be called with the cfg input parameters.
%          Path doesn NOT need to be included.
%
% Example:
%  >> cfg=[]
%  >> cfg.side='r';
%  >> check_cfg(cfg,'elec2parc.m');
%
% Author: David Groppe
% Mehtalab 2013


function checkCfg(cfg,mfile)

if ~isempty(cfg)
    if isempty(mfile),
        error('mfile is empty');
    end
    
    cmnd=sprintf('full_mfile=which(''%s'');',mfile);
    eval(cmnd);
    
    if isempty(full_mfile)
        error('%s not in MATLAB path.',mfile);
    end
        
    % Read in possible fields from mfile
    fid=fopen(full_mfile,'r');
    header=1;
    while header
        % read in line
        new_ln=fgetl(fid);
        
        % remove spaces
        space_ids=find(new_ln==' ');
        n_char=length(new_ln);
        new_ln=new_ln(setdiff(1:n_char,space_ids));
        
        if (length(new_ln)>=15) && strcmpi(new_ln(1:15),'if~isfield(cfg,')
            header=0;
        end
    end
    
    cfg_lines=1;
    mfile_fields=[];
    field_ct=0;
    while cfg_lines
        if (length(new_ln)>=15) && strcmpi(new_ln(1:15),'if~isfield(cfg,')
            quote_ids=find(new_ln==39);
            field_ct=field_ct+1;
            mfile_fields{field_ct}=new_ln((quote_ids(1)+1):quote_ids(2)-1);
        else
            %done reading cfg lines
            cfg_lines=0;
        end
        
        % read in line
        new_ln=fgetl(fid);
        
        % remove spaces
        space_ids=find(new_ln==' ');
        n_char=length(new_ln);
        new_ln=new_ln(setdiff(1:n_char,space_ids));
    end
    fclose(fid);
    
    % Collect any fields of the cfg variable that are not mentioned in the
    % mfile
    cfg_fields=fields(cfg);
    n_cfg_fields=length(cfg_fields);
    found_fields=zeros(1,n_cfg_fields);
    for a=1:n_cfg_fields,
        if ismember(cfg_fields{a},mfile_fields)
            found_fields(a)=1;
        end
    end
    
    non_found_ids=find(found_fields==0);
    if ~isempty(non_found_ids)
        fprintf('Checking cfg fields for function %s\n',full_mfile);
        fprintf('cfg variable has the following fields that are not recognized by the function:\n');
        for a=non_found_ids,
            fprintf('->%s\n',cfg_fields{a});
        end
        error('Remove unrecognized fields from cfg variable.');
    end 
end