function [ccimage] = sbx_ccim_registered(varargin)
% Routine for alignement and 2D-xcoorelation of G/R datasets. Stores
% aligned red, green and xcorr images for ROI selection, together with the
% tforms for later alignement of the imported ROIs.
% 
% Arguments:
% folder = string
% maxshift = [n,n] -> Maximum translation x/y direction
% rotation = [r,r] -> Maximum rotation angle
% align_channel = 'red','green'
% xcorr_channel = 'red','green'

args=struct('folder',[],'maxshift',[50,50],'rotation',[-5:0.5:5],'margin',100,...
    'align_channel','green','xcorr_channel','green','nonrigid',true);

% Overwrite default parameters if required
for pair = reshape(varargin,2,[])
    if isfield(args,pair{1})
        args.(pair{1})=pair{2};
    else
        error('Input argument does not exist\n')
    end
end

if args.folder
    oldfolder = cd(args.folder);
end

if args.nonrigid
    GreenFiles = ls('*_nonrigid.sbx');
else
    GreenFiles = ls('*.align');
end

handles={};

% 1) Create handles for all complete datasets in folder
for fi = size(GreenFiles,1):-1:1
    if args.nonrigid
        filename = strrep(GreenFiles(fi,:),'.sbx','');
    else
        % strip 'green_corr.tif' from filename
        filename = strrep(GreenFiles(fi,:),'.align','');
    end
    
    % Load image
    try
        g = sbximport(filename,1,~args.nonrigid);
        r = sbximport(filename,2,~args.nonrigid);
        
        handles{fi}.green = mean(g,3);
        handles{fi}.red = mean(r,3);
        
        % Create xcorr image
        if strcmp(args.xcorr_channel,'green')
            handles{fi}.ccim=CrossCorrImage(single(g));
        elseif strcmp(args.xcorr_channel,'red')
            handles{fi}.ccim=CrossCorrImage(single(r));
        else
            warning(['2D-Cross correlation failed']);
            return;
        end
        disp([filename ' has been processed.']);
    catch
        warning([filename ' couldnt be processed.']);
    end
end

% 2) Register images, create AVG stacks and tforms
Atform = [1 0 0; 0 1 0; 0 0 1];
tforms{1}=maketform('affine',Atform); % no translation for first image.

greenstack=zeros(size(handles{1}.green,1),size(handles{1}.green,2),length(handles),'uint16');
greenstack(:,:,1)=handles{1}.green;

redstack=zeros(size(handles{1}.red,1),size(handles{1}.red,2),length(handles),'uint16');
redstack(:,:,1)=handles{1}.red; %OMITTED IMCONVERT_UINT8 160625

ccimstack=zeros(size(handles{1}.ccim,1),size(handles{1}.ccim,2),length(handles),'uint16');
ccimstack(:,:,1)=handles{1}.ccim*32768; % Inflate the correlation value (0-1) double to match uint16 range

for m=length(handles):-1:2
    bounds{1} = 1+args.margin:size(handles{1}.green,1)-args.margin;
    bounds{2} = 1+args.margin:size(handles{1}.green,2)-args.margin; % Exclude deadband areas
    
    if strcmp(args.align_channel,'green')
        [~, tforms{m}] = xcorr2_reg(handles{1}.green(bounds{1},bounds{2}),...
            handles{m}.green(bounds{1},bounds{2}),args.maxshift,args.rotation);
    elseif strcmp(args.align_channel,'red')
        [~, tforms{m}] = xcorr2_reg(handles{1}.red(bounds{1},bounds{2}),...
            handles{m}.red(bounds{1},bounds{2}),args.maxshift,args.rotation);
    else
        warning('Specify proper channel for alignment!');
        return
    end
    greenstack(:,:,m)=imtransform(handles{m}.green, tforms{m},...
        'XData',[1 size(handles{1}.green,2)],'YData',[1 size(handles{1}.green,1)]);
    redstack(:,:,m)=imtransform(handles{m}.red, tforms{m},...
        'XData',[1 size(handles{1}.red,2)],'YData',[1 size(handles{1}.red,1)]);
    ccimstack(:,:,m)=imtransform(handles{m}.ccim, tforms{m},...
        'XData',[1 size(handles{1}.ccim,2)],'YData',[1 size(handles{1}.ccim,1)])*32768; % double to uint16
    disp([GreenFiles(m,:) ' has been aligned.']);
end

implay(greenstack);

% 4) Create the average ccimage for ROI identification, save it to file,
% return.
ccimage=mean(ccimstack,3);

SaveWrite_MultiTiff(ccimstack,strcat(args.folder,'MeanCrossCorrImage.tif'));
SaveWrite_MultiTiff(greenstack,strcat(args.folder,'AVGgreenImage.tif'));
SaveWrite_MultiTiff(redstack,strcat(args.folder,'AVGredImage.tif'));
save(strcat(args.folder,'tforms.mat'),'tforms');
end
    
    