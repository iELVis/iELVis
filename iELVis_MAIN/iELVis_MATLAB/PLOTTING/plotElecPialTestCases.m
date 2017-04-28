%%
cfg=[];
cfg.view='r';
cfg.figid=1;
cfg.showlabels='n';
cfg.overlay_parcellation='DK';
cfg.rotate3d='n';
cfg.title='TWH14: Desikan-Killiany Atlas';
cfg_out=plotElecPial('TWH014',cfg);


%%
elecnames=cell(6,1);
for a=1:6,
   elecnames{a}=sprintf('RMF%d',a); 
end
cfg=[];
cfg.view='r';
cfg.figid=1;
cfg.elecshape='sphere';
cfg.eleccolors=rand(6,1);
cfg.colorscale='minmax';
cfg.showlabels='n';
cfg.units='r';
cfg.elecnames=elecnames;
cfg.rotate3d='n';
cfg.title='TWH014: Pie Man Correlations';
cfg_out=plotElecPial('TWH014',cfg);


%%
elecnames=cell(6,1);
for a=1:6,
   elecnames{a}=sprintf('RMF%d',a); 
end
%cfg=[];
cfg.view='lomni';
cfg.figid=2;
cfg.elecshape='sphere';
%cfg.eleccolors=randn(6,1);
cfg.units='t';
cfg.elecnames=elecnames;
cfg.rotate3d='n';
cfg.title=[];
cfg_out=plotElecPial('TWH014',cfg);


%%
cfg=[];
cfg.view='r';
cfg.figid=2;
cfg.rotate3d='n';
cfg.title=[];
cfg.surftype='inflated';
cfg_out=plotElecPial('TWH014',cfg);

%%
cfg=[];
cfg.view='romni';
cfg.figid=3;
cfg.overlay_parcellation='DK';
cfg.rotate3d='n';
cfg.title=[];
cfg_out=plotElecPial('TWH014',cfg);

%%
cfg=[];
cfg.view='romni';
cfg.figid=3;
cfg.surftype='inflated';
cfg.rotate3d='n';
cfg.title=[];
cfg_out=plotElecPial('TWH014',cfg);

%%
cfg=[];
cfg.view='lomni';
cfg.figid=3;
cfg.overlay_parcellation='DK';
cfg.rotate3d='n';
cfg.title=[];
cfg_out=plotElecPial('TWH014',cfg);

%%
cfg=[];
cfg.view='omni';
cfg.figid=3;
%cfg.overlay_parcellation='DK';
cfg.rotate3d='n';
cfg.title=[];
cfg_out=plotElecPial('TWH014',cfg);
