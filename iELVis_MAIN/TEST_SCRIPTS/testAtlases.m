printEm=0;

%%
cfg=[];
cfg.view='l';
cfg.overlayParcellation='DK';
cfg.title='PT001: DK Atlas'; 
cfgOut=plotPialSurf('PT001',cfg);
if printEm
   print(gcf,'-djpeg','pt001atlasDKnoLabels'); 
end

%%
parcOut=elec2Parc('PT001','DK');

%%
cfg=[];
cfg.view='l';
cfg.overlayParcellation='D';
cfg.title='PT001: Destrieux Atlas';
cfgOut=plotPialSurf('PT001',cfg);
if printEm
   print(gcf,'-djpeg','pt001atlasDleftLat'); 
end

%%
parcOut=elec2Parc('PT001','D');

%%
createIndivYeoMapping('PT001');

%%
cfg=[];
cfg.view='l';
cfg.overlayParcellation='Y17';
cfg.title='PT001: Yeo 17-Area';
cfgOut=plotPialSurf('PT001',cfg);
if printEm
   print(gcf,'-djpeg','pt001atlasY17'); 
end


%%
parcOut=elec2Parc('PT001','Y7');


%%
% Read annotation file
[averts, label, col]=read_annotation(fullfile(getFsurfSubDir(),'fsaverage','label','lh.Yeo2011_17Networks_N1000.annot'));

% Network you want to plot
id = 17;

% Make verteces gray
parc_col = .7.*255.*ones(size(col.table(:,1:3)));

% Color parcel of interest
parc_col(id,:)=col.table(id,1:3);

cfg=[];
cfg.view='l';
cfg.overlayParcellation='Y17';
cfg.title=sprintf('Y17 Atlas: Network %d',id); 
cfg.parcellationColors = parc_col;
cfgOut=plotPialSurf('PT001',cfg);
if printEm
   print(gcf,'-djpeg','Yeo17_NW16'); 
end

%%
disp('Script testAtlases.m completed successfully.')