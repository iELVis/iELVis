function plotMgridOnPial(fsub,printEm)
%function plotMgridOnPial(fsub,printEm)
%
% Makes romni or lomni plots on a gray pial and DK atlas pial of a
% patient's brain with electrodes colored according to mgrid file colors.
%
% Inputs:
%  fsub  - FreeSurfer subject directory
%  printEm - If nonzero jpgs of the figures will be made in the subject's
%            elec_recon folder
%
% Examples:
% plotMgridOnPial('TWH10',0);
%
%
% Author:
% David M. Groppe
% March, 2015
%

% Future Work:
%  -Add an option to make depths invisible

% History:
% May 2015-Now electrode pairs are derived directly from mgrid file and
% disabled electrodes are not shown.

% FreeSurfer Subject Directory
fsDir=getFsurfSubDir();

if nargin<2,
    printEm=0;
end


%% Get mgrid info
%[~, elecLabels, elecRgb]=mgrid2matlab(fsub,hem);
[~, elecLabels, elecRgb, elecPairs, elecPresent]=mgrid2matlab(fsub);
nElec=length(elecLabels);

elecnames=cell(1,length(elecLabels));
elecHem=zeros(nElec,1); %1=left
elecDepth=zeros(nElec,1);
elecType=cell(nElec,1);

for a=1:length(elecLabels),
    elecnames{a}=rmChar(elecLabels{a}(3:end),'_');
    if strcmpi(elecLabels{a}(1),'L')
        elecHem(a)=1;
    end
    elecType{a}=elecLabels{a}(2);
%     if strcmpi(elecType{a},'D')
%         elecDepth(a)=1;
%     end
end

pairPresent=zeros(size(elecPairs,1),1);
for a=1:size(elecPairs,1),
    elecPairs{a,4}=lower(elecPairs{a,1}(1));
    elecPairs{a,1}=rmChar(elecPairs{a,1}(3:end),'_');
    elecPairs{a,2}=rmChar(elecPairs{a,2}(3:end),'_');
    elecId1=findStrInCell(elecPairs{a,1},elecnames,1);
    elecId2=findStrInCell(elecPairs{a,2},elecnames,1);
    pairPresent(a)=elecPresent(elecId1)*elecPresent(elecId2);
end

%% Get Unique Electrode Names and Colors for Legend
uniStemsL=[];
uniStemsRgbL=[];
leftIds=find(elecHem);
ct=0;
for a=leftIds',
    underIds=find(elecLabels{a}=='_');
    elecstem=elecLabels{a}(underIds(1)+1:underIds(2)-1);
    
    if ~ismember(elecstem,uniStemsL)
        ct=ct+1;
        uniStemsL{ct}=elecstem;
        uniStemsRgbL(ct,1:3)=elecRgb(a,:);
    end 
end

uniStemsR=[];
uniStemsRgbR=[];
rightIds=find(~elecHem);
ct=0;
for a=rightIds',
    underIds=find(elecLabels{a}=='_');
    elecstem=elecLabels{a}(underIds(1)+1:underIds(2)-1);
    
    if ~ismember(elecstem,uniStemsR)
        ct=ct+1;
        uniStemsR{ct}=elecstem;
        uniStemsRgbR(ct,1:3)=elecRgb(a,:);
    end 
end


%%
elecPresentIds=find(elecPresent);
for hemLoop=1:2,
    if hemLoop==1,
        hem='l';
        hemLong='Left';
        hemIds=leftIds;
        uniStems=uniStemsL;
        uniStemsRgb=uniStemsRgbL;
    else
        hem='r';
        hemLong='Right';
        hemIds=rightIds;
        uniStems=uniStemsR;
        uniStemsRgb=uniStemsRgbR;
    end
    useIds=intersect(elecPresentIds,hemIds);

    if ~isempty(hemIds)
        for fLoop=1:2,
            cfg=[];
            cfg.view=[hem 'omni'];
            %cfg.figId=fLoop+(hemLoop-1)*10;
            cfg.elecColors=elecRgb(useIds,:);
            cfg.elecNames=elecnames(useIds);
            cfg.elecCbar='n';
            cfg.ignoreDepthElec='n';
            cfg.title=fsub;
            %cfg.title=[];
            cfg.pairs=elecPairs;
            %cfg.pairs=elecPairs(useIds,:);
            %cfg.showlabels='y';
            if fLoop==2
                cfg.overlayParcellation='DK';
            end
            %cfg.rotate3d='n';
            cfgOut=plotPialSurf(fsub,cfg);
            
            % Prevent the code below from creating a legend on the previous
            % plot
            pause(.5);
            drawnow;
            
            % Add electrode legend
            hAx=axes('position',[.9 .01 .07 .98]);
            v=axis;
            nUni=length(uniStems);
            dlt=(v(4)-v(3))/(nUni+1);
            for a=1:nUni,
                ht=text(v(2),v(3)+dlt*a,uniStems{a});
                set(ht,'fontsize',18,'color',uniStemsRgb(a,:), ...
                    'horizontalalignment','right','fontweight','bold', ...
                    'backgroundcolor','k');
            end
            set(hAx,'box','off','visible','off');
            
            set(gcf,'paperpositionmode','auto');
            if universalYes(printEm)
                % Make sure PICS directory exists
                erPath=fullfile(fsDir,fsub,'elec_recon');
                outPath=fullfile(erPath,'PICS');
                if ~exist(outPath,'dir')
                    dirSuccess=mkdir(outPath);
                    if ~dirSuccess,
                        error('Could not create directory %s',dirSuccess);
                    end
                end
                
                if fLoop==1,
                    figFname=sprintf('%s%sMgridElec',fsub,hemLong);
                else
                    figFname=sprintf('%s%sMgridElecDK',fsub,hemLong);
                end
                outFname=fullfile(erPath,'PICS',figFname);
                print(gcf,outFname,'-djpeg');
            end
        end
    end
end
