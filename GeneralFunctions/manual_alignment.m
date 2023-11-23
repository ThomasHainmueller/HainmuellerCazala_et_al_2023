function manual_alignment(varargin)
% Use a folder with preexisting tforms and go through recordings/channels
% one by one to manually correct alignment errors.

global flag
flag = 1;

args=struct('alignment_channel','green','selection','all','rois','saved'); 
% planes: e.g [1,2,3]; channels: 'green','red','ccim'; rois: 'saved',
% 'previous'

% Overwrite default parameters if required
for pair = reshape(varargin,2,[])
    if isfield(args,pair{1})
        args.(pair{1})=pair{2};
    else
        error('Input argument does not exist\n')
    end
end

switch args.alignment_channel
    case 'green'
        chan = 1;
    case 'red'
        chan = 2;
    case 'ccim'
        chan = 3;
end

fnames = dir('*.sbx');

if strcmp(args.selection,'all')
    args.selection = 1:length(fnames);
end

Ref = load(strrep(fnames(1).name,'.sbx','.mat'),'images');
Ref = Ref.images;

for index = 1:length(args.selection)
    fi = args.selection(index);
    fn = strrep(fnames(fi).name,'.sbx','.mat');

    % Load initial tform, previous is usefull for coplex deformations
    if strcmp(args.rois,'saved') || index==1
        tform = load(fn,'tform');
    elseif strcmp(args.rois,'previous')
        tn = strrep(fnames(args.selection(index-1)).name,'.sbx','.mat');
        tform = load(tn,'tform');
    end
    
    tform = tform.tform;
        
    Mov = load(fn,'images');
    Mov = Mov.images;
    
    for pl = 1:size(Ref,1)
        if fi == 1
            for ch = 1:size(Ref,4)
                imwrite(uint16(squeeze(Ref(pl,:,:,ch))),...
                    sprintf('Mean_plane%d_channel%d.tif',pl,ch),...
                    'TIFF','writemode','append')
            end
        else
            flag = 0; % Make sure windows are called one after the other
            
            % Core of algorithm, correct the tform
            aligndisplay(squeeze(Ref(pl,:,:,chan)),...
                squeeze(Mov(pl,:,:,chan)), fn, pl, tform{pl});
            
            while ~flag
                pause(.1); % Wait till GUI returns
            end
            
            % transform image, append to tiff stack
            for ch = 1:size(Ref,4)
                % Reload tform after modification
                tform = load(fn,'tform');
                tform = tform.tform;
                
                thisMov = imtransform(squeeze(Mov(pl,:,:,ch)),tform{pl},...
                    'XData',[1 size(Ref,3)], 'YData', [1 size(Ref,2)]); %CAVE: Check whether X/Y is correct here!
                
                imwrite(uint16(thisMov),...
                    sprintf('Mean_plane%d_channel%d.tif',pl,ch),...
                    'TIFF','writemode','append')
            end            
        end
    end
end

function [] = aligndisplay(A,B,filename,plane,tform)
% Display two images via imfuse. Allow the user to shift and rotate B
% relative to A to achieve alignment. Save modified tform to specified file
% 
% CAVE: This code is highly specific due to the limited utility of the
% figure handle/arrow buttons. Therefore it only works with a complete
% dataset processed with motion corr in the standard manner 170808

global a Aim Bim fn pl
Aim = A;
fn = filename;

% Iniate with no translation if no tform is provided
if nargin < 5
    a = [1,0,0; 0,1,0; 0,0,1];
    tform = maketform('affine',a);
else
    a = tform.tdata.T;
end

pl = plane;
%Bim = imtransform(B,tform,'XData',[1 size(Aim,1)],'YData',[1 size(Aim,2)]);
Bim = B;

S.fh = figure('keypressfcn',@fh_kpfcn,'Position',[10,10,1900,1200]);
S.tx = uicontrol('style','text',...
                  'units','pixels',...
                  'position',[60 120 80 20],...
                  'fontweight','bold'); 

key = 0;
guidata(S.fh,S)


function [] = fh_kpfcn(H,E)
global a Aim Bim fn pl flag

% Figure keypressfcn
S = guidata(H);
P = get(S.fh,'position');
set(S.tx,'string',E.Key)
switch E.Key
    
    % translation
    case {'e','uparrow'}
        a(3,2) = a(3,2)-1;
    case {'d','downarrow'}
        a(3,2) = a(3,2)+1;
    case {'s','leftarrow'}
        a(3,1) = a(3,1)-1;
    case {'f','rightarrow'}
        a(3,1) = a(3,1)+1;
    
    % rotation
    case 'w'
        m = [cosd(.2) sind(.2) 0; -sind(.2) cosd(.2) 0; 0 0 1];
        a = a*m;
%         r = asind(a(1,2))+.2;
%         sc = a(1,1) - cosd(asind(a(1,2))) + cosd(r);
%         a = [sc sind(r) 0; -sind(r) sc 0; a(3,1) a(3,2) 1];
%         %a = [cosd(r) sind(r) 0; -sind(r) cosd(r) 0; a(3,1) a(3,2) 1];
    case 'r'
        m = [cosd(.2) -sind(.2) 0; sind(.2) cosd(.2) 0; 0 0 1];
        a = a*m;
%        r = asind(a(1,2))-.2;
%        sc = a(1,1) - cosd(asind(a(1,2))) + cosd(r);
%        a = [sc sind(r) 0; -sind(r) sc 0; a(3,1) a(3,2) 1];       
%        %a = [cosd(r) sind(r) 0; -sind(r) cosd(r) 0; a(3,1) a(3,2) 1];
    
    % scaling
    case 'l'
        a(1,1) = a(1,1)+.005;
    case 'j'
        a(1,1) = a(1,1)-.005;
    case 'i'
        a(2,2) = a(2,2)+.005;
    case 'k'
        a(2,2) = a(2,2)-.005;
    
    %return
    case {'b', 'return'}
        % Cave: image must be moved at least once before exiting
        tform = load(fn,'tform');
        tform = tform.tform;
        tform{pl} = maketform('affine',a);
        save(fn,'tform','-append');
        flag = 1;
        close all
        fprintf('%s, plane %d saved\n',fn,pl);
        return
end

tform = maketform('affine',a);
thisBim = imtransform(Bim,tform,'XData',[1 size(Aim,2)],'YData',[1 size(Aim,1)]);
disp = imfuse(Aim,thisBim);
image(disp);
drawnow;
