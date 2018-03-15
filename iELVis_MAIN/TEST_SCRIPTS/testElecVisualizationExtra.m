
%% Plot electrodes color coded to represent correlations OMNI view, Unilateral coverage
elecNames=cell(8,1);
for a=1:8,
    elecNames{a}=sprintf('LGd%d',a+16);
end
cfg=[];
cfg.view='omni';
cfg.figId=1;
% cfg.elecShape='sphere';
cfg.elecColors=rand(8,1);
cfg.elecColorScale='minmax';
cfg.showLabels='n';
cfg.elecUnits='r';
cfg.elecNames=elecNames;
cfg.title='PT001: Stimulus Correlations';
cfgOut=plotPialSurf('PT001',cfg);


%% Plot electrodes color coded to represent correlations LOMNI view, Unilateral coverage
elecNames=cell(8,1);
for a=1:8,
    elecNames{a}=sprintf('LGd%d',a+16);
end
cfg=[];
cfg.view='lomni';
cfg.figId=1;
% cfg.elecShape='sphere';
%cfg.elecColors=rand(8,1);
cfg.elecColors=linspace(0,1,8)';
cfg.elecColorScale='minmax';
cfg.showLabels='n';
cfg.elecUnits='r';
cfg.elecNames=elecNames;
cfg.title='PT001: Stimulus Correlations';
cfgOut=plotPialSurf('PT001',cfg);


%% Plot electrodes color coded to represent correlations ROMNI view, Unilateral coverage on the 
% non-visualized hemisphere
% elecNames=cell(8,1);
% for a=1:8,
%     elecNames{a}=sprintf('LGd%d',a+16);
% end
% cfg=[];
% cfg.view='romni';
% % cfg.view='r';
% cfg.figId=1;
% % cfg.elecShape='sphere';
% %cfg.elecColors=rand(8,1);
% cfg.elecColors=linspace(0,1,8)';
% cfg.elecColorScale='minmax';
% cfg.showLabels='n';
% cfg.elecUnits='r';
% cfg.elecNames=elecNames;
% cfg.title='PT001: Stimulus Correlations';
% cfgOut=plotPialSurf('PT001',cfg);
%
% This now throws an error as it should

%% Inlfated Brain, Lomni view
cfg=[];
cfg.view='lomni';
cfg.figId=1;
cfg.surfType='inflated';
cfg.title=[];
cfgOut=plotPialSurf('PT001',cfg);


%% Plot electrodes color coded to represent correlations OMNI view, Bilateral coverage
sub='TWH013';
elecNames=cell(12,1);
ct=0;
for a=1:8,
    ct=ct+1;
    elecNames{ct}=sprintf('LTGRID%d',a+16);
end
for a=1:4,
    ct=ct+1;
    elecNames{ct}=sprintf('RPT%d',a);
end
cfg=[];
cfg.view='omni';
cfg.figId=1;
% cfg.elecShape='sphere';
cfg.elecColors=rand(ct,1);
cfg.elecColorScale='minmax';
cfg.showLabels='n';
cfg.elecUnits='r';
cfg.elecNames=elecNames;
cfg.elecSize=6;
cfg.title=sprintf('%s: Stimulus Correlations',sub);
cfgOut=plotPialSurf(sub,cfg);


%% Plot electrodes color coded to represent correlations LOMNI view, Bilateral coverage
sub='TWH013';
elecNames=cell(12,1);
ct=0;
for a=1:8,
    ct=ct+1;
    elecNames{ct}=sprintf('LTGRID%d',a+16);
end
for a=1:4,
    ct=ct+1;
    elecNames{ct}=sprintf('RPT%d',a);
end
cfg=[];
cfg.view='lomni';
cfg.figId=1;
% cfg.elecShape='sphere';
cfg.elecColors=rand(ct,1);
cfg.elecColorScale='minmax';
cfg.showLabels='n';
cfg.elecUnits='r';
cfg.elecNames=elecNames;
cfg.title=sprintf('%s: Stimulus Correlations',sub);
cfgOut=plotPialSurf(sub,cfg);


%% Plot electrodes color coded to represent correlations ROMNI view, Bilateral coverage
sub='TWH013';
elecNames=cell(12,1);
ct=0;
for a=1:8,
    ct=ct+1;
    elecNames{ct}=sprintf('LTGRID%d',a+16);
end
for a=1:4,
    ct=ct+1;
    elecNames{ct}=sprintf('RPT%d',a);
end
cfg=[];
cfg.view='lomni';
cfg.figId=1;
% cfg.elecShape='sphere';
cfg.elecColors=rand(ct,1);
cfg.elecColorScale='minmax';
cfg.showLabels='n';
cfg.elecUnits='r';
cfg.elecNames=elecNames;
cfg.title=sprintf('%s: Stimulus Correlations',sub);
cfgOut=plotPialSurf(sub,cfg);



%% Plot electrodes color coded to represent correlations Left Lat view, Bilateral coverage
sub='TWH013';
elecNames=cell(12,1);
ct=0;
for a=1:8,
    ct=ct+1;
    elecNames{ct}=sprintf('LTGRID%d',a+16);
end
for a=1:4,
    ct=ct+1;
    elecNames{ct}=sprintf('RPT%d',a);
end
cfg=[];
cfg.view='l';
cfg.figId=1;
% cfg.elecShape='sphere';
cfg.elecColors=rand(ct,1);
cfg.elecColorScale='minmax';
cfg.showLabels='n';
cfg.elecUnits='r';
cfg.elecNames=elecNames;
cfg.title=sprintf('%s: Stimulus Correlations',sub);
cfgOut=plotPialSurf(sub,cfg);

disp('Ignore the warning message just thrown by plotPialSurf. Right hem electrodes were included as a bug check in this example.');
disp('It is supposed to trigger a warning.');


%% Plot a single electrode with color specified (this produced an error for Kathrin once)
elecNames=cell(1,1);
for a=1:1,
    elecNames{a}=sprintf('LGd%d',a+16+4);
end
cfg=[];
cfg.view='l';
cfg.elecShape='sphere';
cfg.elecColors=[1 0 1];
cfg.elecColorScale=[.5 .5];
cfg.elecCbar='n';
cfg.elecSize=3;
cfg.showLabels='n';
cfg.elecUnits='r';
cfg.elecNames=elecNames;
cfg.title='PT001: Stimulus Correlations';
cfgOut=plotPialSurf('PT001',cfg);


%%
disp('testElecVisualizationExtra.m completed successfully.');


return;

%%%%%%% ERROR CASES %%%%%%%

%% Throw an error when multiple electrodes have the exact same name
elecNames=cell(8,1);
for a=1:8,
    elecNames{a}=sprintf('LGd1');
end
cfg=[];
cfg.view='lomni';
cfg.elecColors=rand(8,1);
cfg.elecColorScale='minmax';
cfg.showLabels='n';
cfg.elecUnits='r';
cfg.elecNames=elecNames;
cfg.title='PT001: Stimulus Correlations';
cfgOut=plotPialSurf('PT001',cfg);

