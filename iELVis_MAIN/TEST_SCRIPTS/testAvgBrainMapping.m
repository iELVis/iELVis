%%
printEm=0;

%% Test Wiki Example 1
groupAvgCoords=[];
groupLabels=[];
groupIsLeft=[];
subs={'PT001','PT002'};
cfg=[];
% cfg.rmDepths=1;
for a=1:length(subs),
    fprintf('Working on Participant %s\n',subs{a});
    [avgCoords, elecNames, isLeft]=sub2AvgBrain(subs{a},cfg);
    groupAvgCoords=[groupAvgCoords; avgCoords];
    groupLabels=[groupLabels; elecNames];
    groupIsLeft=[groupIsLeft; isLeft];
end

%% Test Wiki Example 2
cfg=[];
cfg.view='l';
cfg.elecCoord=[groupAvgCoords groupIsLeft];
cfg.elecNames=groupLabels;
cfg.showLabels='n';
cfg.title='PT001-2 on Avg. Brain';
cfgOut=plotPialSurf('fsaverage',cfg);
if printEm,
    print(gcf,'-djpeg','pt001-2onAvgBrain');
end

%% Same as above but on Inflated average brain
cfg=[];
cfg.view='l';
cfg.elecCoord=[groupAvgCoords groupIsLeft];
cfg.elecNames=groupLabels;
cfg.showLabels='n';
cfg.surfType='inflated';
cfg.title='PT001-2 on Avg. Brain';
cfgOut=plotPialSurf('fsaverage',cfg);


%% Test MNI305 mapping
[avgCoords, elecNames]=subElecs2MNI305('PT001');


%%
disp('Script testAvgBrainMapping.m completed successfully.')