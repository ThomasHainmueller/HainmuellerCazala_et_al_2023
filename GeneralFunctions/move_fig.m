function [] = move_fig()
% move figure with arrow keys.

S.fh = figure('keypressfcn',@fh_kpfcn);
% S.tx = uicontrol('style','text',...
%                  'units','pixels',...
%                  'position',[60 120 80 20],...
%                  'fontweight','bold'); 
guidata(S.fh,S)             

function [a] = fh_kpfcn(H,E)          
% Figure keypressfcn
S = guidata(H);
P = get(S.fh,'position');
set(S.tx,'string',E.Key)
switch E.Key
    case 'rightarrow'
        set(S.fh,'pos',P+[5 0 0 0])
        a=1
    case 'leftarrow'
        set(S.fh,'pos',P+[-5 0 0 0])
        a=2
    case 'uparrow'
        set(S.fh,'pos',P+[0 5 0 0])
    case 'downarrow'
        set(S.fh,'pos',P+[0 -5 0 0])
    case 'return'
        close(S.fh);
        return
    otherwise  
end
