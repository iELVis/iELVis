%% Analogous to Test Wiki Example 1, but with Depths
groupAvgCoords=[];
groupLabels=[];
groupIsLeft=[];
subs={'NiAs','AnRo'};
cfg=[];
% cfg.rmDepths=1;
for a=1:length(subs),
    fprintf('Working on Participant %s\n',subs{a});
    [avgCoords, elecNames, isLeft]=sub2AvgBrain(subs{a},cfg);
    groupAvgCoords=[groupAvgCoords; avgCoords];
    groupLabels=[groupLabels; elecNames];
    groupIsLeft=[groupIsLeft; isLeft];
end

%% Analogous to Test Wiki Example 2, but with Depths
cfg=[];
cfg.view='r';
cfg.elecCoord=[groupAvgCoords groupIsLeft];
cfg.elecNames=groupLabels;
cfg.showLabels='n';
cfg.opaqueness=0.5;
cfg.title='NiAs & AnRo on Avg. Brain';
cfgOut=plotPialSurf('fsaverage',cfg);


%% Same as above but on Inflated average brain
cfg=[];
cfg.view='r';
cfg.elecCoord=[groupAvgCoords groupIsLeft];
cfg.elecNames=groupLabels;
cfg.showLabels='n';
cfg.surfType='inflated';
cfg.title='PT001-2 on Avg. Brain';
cfgOut=plotPialSurf('fsaverage',cfg);


%%
disp('Script testAvgBrainMappingExtra.m completed successfully.')