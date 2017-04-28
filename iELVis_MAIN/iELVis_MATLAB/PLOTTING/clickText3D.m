function clickText3D(h,txt,pop_fact)
%function clickText3D(h,txt,pop_fact)
%
% Same as click_text.m, but when you click on the object, the text appears
% out towards the viewer so that it doesn't get obscured by other objects
% (e.g., the brain) in the plot.
%
% Inputs:
%  h        - The handle of an object in a figure.
%  txt      - The text that will be displayed when you click on that object.
%  pop_fact - The factor by which the text will be projected out towards
%             the viewer. 0 means the text will appear on the object. The 
%             more positive the value, the further the text will be 
%             projected towards the viewer.
%
% Author:
% David Groppe

set(h,'userdata',txt);
bdfcn=['Cp = get(gca,''CurrentPoint''); ' ...
    'Cp=Cp(1,1:3);', ...
    'v=axis;', ...
    'campos=get(gca,''cameraposition'');', ...
    'df=Cp-campos;', ...
    'nrmd=df/sqrt(sum(df.^2));', ...
    sprintf('Cp=Cp-%d*nrmd;',pop_fact), ...
    'dat=get(gcbo,''userdata'');', ...
    'ht=text(Cp(1),Cp(2),Cp(3),sprintf(''%s'',dat));', ...
    'set(ht,''backgroundcolor'',''w'',''horizontalalignment'',''center'',''verticalalignment'',''middle'',''buttondownfcn'',''delete(gcbo);'');'];
set(h,'buttondownfcn',bdfcn);
