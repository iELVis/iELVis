%% Monopolar
elecnames=[];
for a=1:4,
   elecnames{a}=sprintf('LOF%d',a);
end
for a=1:4,
   elecnames{length(elecnames)+1}=sprintf('LLF%d',a);
end
for a=1:4,
   elecnames{length(elecnames)+1}=sprintf('LAT%d',a);
end
cfg=[];
cfg.view='l';
cfg.figid=1;
cfg.elecshape='sphere';
%cfg.elecshape='marker';
cfg.eleccolors=rand(length(elecnames),1)*2-1;
cfg.colorscale='minmax';
cfg.showlabels='n';
cfg.units='r';
cfg.elecnames=elecnames;
cfg.rotate3d='n';
cfg.title='TiFl: Dummy Data';
cfg_out=plotElecPial('TiFl',cfg);


%% Bipolar
pairs=[];
ct=0;
for a=1:3,
    ct=ct+1;
    pairs{ct,1}=sprintf('LPT%d',a);
    pairs{ct,2}=sprintf('LPT%d',a+1);
    pairs{ct,3}=rand(1,3); % RGB val
    pairs{ct,4}='L';
end
for a=1:3,
    ct=ct+1;
    pairs{ct,1}=sprintf('LLF%d',a);
    pairs{ct,2}=sprintf('LLF%d',a+1);
    pairs{ct,3}=rand(1,3); % RGB val
    pairs{ct,4}='L';
end
cfg=[];
cfg.view='l';
cfg.figid=2;
cfg.elecshape='sphere';
%cfg.elecshape='marker';
cfg.pairs=pairs;
%cfg.eleccolors=rand(length(elecnames),1)*2-1;
cfg.colorscale='minmax';
cfg.showlabels='n';
cfg.units='r';
cfg.elecnames=elecnames;
cfg.rotate3d='n';
cfg.title='TiFl: Dummy Data';
cfg_out=plotElecPial('TiFl',cfg);