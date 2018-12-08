function [figH, axH]=plotPostImpVsLepto(sub,printEm,plotPial)
%function [figH, axH]=plotPostImpVsLepto(sub,printEm,plotPial)
%
% This function creates two plots to illustrate the effect of brain shift
% correction:
%  1: A plot of the distance between each electrode's pre and post
%  correction positions and a 3D scatter plot of these positions
%
%  2: Post-correction positions overlayed on pial surface and color coded
%  to represent distance between pre and post correction locations
%
% This function is called by interpStripElec.m, dykstraElecPjct.m, and yangWangElecPjct.m
%
% Inputs:
%  sub - Subject's FreeSurfer name/folder
%  printEm  - If nonzero, copies of the images are output the
%             "elec_recon/PICS" subfolder of the subject's FreeSurfer folder
%  plotPial - If nonzero, a figure is created that illustrateds
%             post-brainshift corrected electrodes color coded by the extent of
%             brainshift correction.
%
% Author: David Groppe
% Honeylab, Univ. of Toronto
% June 2015

% Load lepto and CT coordinates
fsDir=getFsurfSubDir();
%erPath=[fsDir sub '/elec_recon/'];
erPath=fullfile(fsDir,sub,'elec_recon');
leptoFname=fullfile(erPath,[sub '.LEPTO']);
%leptoFname=[erPath sub '.LEPTO'];
leptoCsv=csv2Cell(leptoFname,' ',2);
nChan=size(leptoCsv,1);

% Load elec names etc..
chanFname=fullfile(erPath,[sub '.electrodeNames']);
chanInfo=csv2Cell(chanFname,' ',2);
chanName=chanInfo(:,1);
chanType=chanInfo(:,2);
chanHem=chanInfo(:,3);

postImpFname=fullfile(erPath,[sub '.POSTIMPLANT']);
postImpCsv=csv2Cell(postImpFname,' ',2);

leptoRAS=zeros(nChan,3);
ctRAS=zeros(nChan,3);
for a=1:nChan,
    for b=1:3,
        ctRAS(a,b)=str2num(postImpCsv{a,b});
        leptoRAS(a,b)=str2num(leptoCsv{a,b});
    end
end


%% Plot results to double check
shiftDist=sqrt( sum( (leptoRAS-ctRAS).^2,2)); %units are mm

%rgb=zeros(nElec,3);
rgb=vals2Colormap(shiftDist,'justpos');
figH(1)=figure;
%set(figH(1),'position',[104 285 1114 410]);
axH(1)=subplot(121);
plot(shiftDist,'.-'); hold on;
last_type=[];
marker='o';
non_depth=ones(1,nChan);
leftHem=ones(1,nChan);
for a=1:nChan,
    if chanType{a}=='D'
        non_depth(a)=0;
        marker='s';
    else
        marker='o';
    end
    if chanHem{a}=='R'
        leftHem(a)=0;
    end
    h=plot(a,shiftDist(a),marker);
    set(h,'color',rgb(a,:));
    clickText(h,[num2str(a) ': ' chanName{a}]);
    xlabel('Channel');
    ylabel('Pre vs. Post Brain Shift Corrected Distance (mm)');
end
non_depth_ids=find(non_depth);
title(sprintf('%s Median=%.1f, SIQR=%.1f (depths ignored)',sub,median(shiftDist(non_depth_ids)), ...
    iqr(shiftDist(non_depth_ids))/2));
v=axis;
axis([1 length(shiftDist) v(3:4)]);

%3D plot of and pre vs. post shift correction locations
axH(2)=subplot(122);
for a=1:nChan,
    h=plot3(leptoRAS(a,1),leptoRAS(a,2),leptoRAS(a,3),'r.'); hold on;
    set(h,'color',rgb(a,:),'markersize',25);
    clickText(h,chanName{a});
    h=plot3(ctRAS(a,1),ctRAS(a,2),ctRAS(a,3),'bo');
    %clickText(h,rm_substring(labels{a},'_'));
    clickText(h,chanName{a});
    plot3([leptoRAS(a,1) ctRAS(a,1)],[leptoRAS(a,2) ctRAS(a,2)],[leptoRAS(a,3) ctRAS(a,3)],'k-');
end
if mean(leftHem)>=0.5,
    % Left electrode majority
    view([-76 12]);
else
    view([100 16]);
end
axis tight;
axis square;
title('Red=Postcorrection, Blue=Precorrection');
xlabel('Left- Right+');
ylabel('Pos- Ant+');
zlabel('Inf- Sup+');

if universalYes(printEm)
    % Make sure PICS directory exists
    outPath=fullfile(erPath,'PICS');
    if ~exist(outPath,'dir')
        dirSuccess=mkdir(outPath);
        if ~dirSuccess,
            error('Could not create directory %s',dirSuccess);
        end
    end
    outFigFname=fullfile(outPath,sprintf('%s_ShiftDist.jpg',sub));
    print(figH(1),'-djpeg',outFigFname);
    outFigFname=fullfile(outPath,sprintf('%s_ShiftDist',sub));
    saveas(figH(1),[outFigFname '.fig']); % PM edit 20160405
end

% Plot shift distances on pial surface
if ~isempty(non_depth_ids)
    if universalYes(plotPial)
        if sum(leftHem)
            if sum(~leftHem)
                vw='omni';
            else
                vw='lomni';
            end
        else
            vw='romni';
        end
        figH(2)=figure;
        cfg=[];
        cfg.view=vw;
        cfg.figId=figH(2);
        cfg.elecColors=shiftDist;
        cfg.elecColorScale='justpos';
        cfg.elecUnits='mm';
        cfg.elecNames=chanName;
        cfg.showLabels='n';
        cfg.title=sprintf('%s: CT to Leptomeningeal distance',sub);
        plotPialSurf(sub,cfg);
        
        if universalYes(printEm)
            outFigFname=fullfile(outPath,sprintf('%s_ShiftDistOnBrain.jpg',sub));
            print(figH(2),'-djpeg',outFigFname);
        end
    end 
end
