printEm=0;

%% DK atlas with elecs and labels
cfg=[];
cfg.view='l';
cfg.overlayParcellation='DK';
cfg.showLabels='y';
cfg.title='PT001: DK Atlas'; 
cfgOut=plotPialSurf('PT001',cfg);
if printEm,
    print(gcf,'-djpeg','pt001atlasDK');
end

%% D atlas with elecs and labels
cfg=[];
cfg.view='omni';
cfg.overlayParcellation='D';
cfg.showLabels='y';
cfg.title='PT001: Destrieux Atlas';
cfgOut=plotPialSurf('PT001',cfg);
if printEm,
    print(gcf,'-djpeg','pt001atlasD');
end 

%% Plot depths with transparent pial surface
cfg=[];
cfg.view='li';
cfg.ignoreDepthElec='n';
cfg.opaqueness=0.5;
cfg.onlyShow={'LDAm1','LDAm2','LDAm3','LDAm4','LDAm5','LDAm6','LDAm7','LDAm8','LDHp1','LDHp2','LDHp3','LDHp4','LDHp5','LDHp6','LDHp7','LDHp8'};
cfg.title='PT001';
cfgOut=plotPialSurf('PT001',cfg);
if printEm,
    print(gcf,'-djpeg','pt001transparent');
end 

%% Plot electrodes color coded to represent correlations
elecNames=cell(8,1);
for a=1:8,
    elecNames{a}=sprintf('LGd%d',a+16);
end
cfg=[];
cfg.view='l';
cfg.elecShape='sphere';
cfg.elecColors=rand(8,1);
cfg.elecColorScale='minmax';
cfg.showLabels='n';
cfg.elecUnits='r';
cfg.elecNames=elecNames;
cfg.elecSize=2;
cfg.title='PT001: Stimulus Correlations';
cfgOut=plotPialSurf('PT001',cfg);
if printEm,
    print(gcf,'-djpeg','pt001corr');
end 


%% Plot electrodes with Bars
pairs=[];
ct=0;
for a=1:7,
  ct=ct+1;
  pairs{ct,1}=sprintf('LGd%d',a+8);
  pairs{ct,2}=sprintf('LGd%d',a+1+8);
  pairs{ct,3}=rand(1,3); % RGB val
  pairs{ct,4}='L';
end
cfg=[];
cfg.view='l';
cfg.pairs=pairs;
cfg.showLabels='n';
cfg.elecUnits='r';
cfg.title='PT001: Stimulus Correlations';
cfg_out=plotPialSurf('PT001',cfg);
if printEm,
    print(gcf,'-djpeg','pt001bipolar');
end 
 
%%
disp('Script testElecVisualization.m completed successfully.')
