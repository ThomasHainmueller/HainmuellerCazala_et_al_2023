function [ccimage] = ccim_registered2(folder, maxshift, rotation, align_channel, xcorr_channel,nodeadbands)
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

if nargin < 6
    nodeadbands=false;
end
if nargin < 5
    xcorr_channel='green'; % Special option for jRGECO
end
if nargin < 4
    align_channel='red';
end
if nargin < 3
    rotation = [0];
end
if nargin < 2
    maxshift=[20,20];
end
if nargin < 1
    folder=strcat(cd(),'\');
end

oldfolder = cd(folder);
GreenFiles = ls('*green_corr.tif');
cd(oldfolder);
handles={};

% 1) Create handles for all complete datasets in folder
for n = size(GreenFiles,1):-1:1
    % strip 'green_corr.tif' from filename
    filename = strrep(GreenFiles(n,:),'green_corr.tif','');
    GreenName = GreenFiles(n,:);
    RedName = strrep(GreenFiles(n,:),'green_corr.tif','red_corr.tif');
    handles{n}.filename=filename;
    
    % Load green and red mean images
    try
        handles{n}.red = mean(ImportMultiTiff(strcat(folder,RedName)),3);
        handles{n}.green = mean(ImportMultiTiff(strcat(folder,GreenName)),3);
    catch
        warning([filename ' couldnt be loaded, incomplete dataset']);
        continue
    end
    
    % Prepare images for alignment
    if strcmp(align_channel,'green')
        %handles{n}.align=handles{n}.green-imopen(handles{n}.green,strel('disk',4));
        handles{n}.align=handles{n}.green;
    elseif strcmp(align_channel,'red')
        handles{n}.align=BGsubt(handles{n}.red,0.001);
    else
        warning(['Alignment channel does not exist']);
        return;
    end
    
    xcorrname = strrep(strcat(folder,GreenFiles(n,:)),'green_corr.tif','ccim.mat');
    % Prepare xcorr images and store mean with the files
    if ~isempty(dir(xcorrname))
        handles{n}.ccim=double(importdata(xcorrname)); % If a ccim was allready generated
    elseif strcmp(xcorr_channel,'green')
        handles{n}.ccim=CrossCorrImage(double(ImportMultiTiff(strcat(folder,GreenName))));
        %handles{n}.ccim=handles{n}.green;
        ccim=handles{n}.ccim;
        save(xcorrname,'ccim');
        clear ccim;
    elseif strcmp(xcorr_channel,'red')
        handles{n}.ccim=CrossCorrImage(double(ImportMultiTiff(strcat(folder,RedName))));
        save(strcat(folder,xcorrname),'handles{n}.ccim');
    else
        warning(['2D-Cross correlation failed']);
        return;
    end
    
    fprintf('%s has been processed.\n',filename);
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
    fprintf('%s is beeing aligned.\n',handles{m}.filename);
    if nodeadbands
        fixed=handles{1}.align;
        moving=handles{m}.align;
        
        ghist = mean(moving,1); % Exclude mirror flyback bands
        thresh = prctile(ghist,1)+0.1*(max(ghist)-prctile(ghist,1));
        b = [find(ghist>thresh,1,'first'),find(ghist>thresh,1,'last')];
        
        [~,tforms{m}] = xcorr2_reg(fixed(:,b(1):b(2)),moving(:,b(1):b(2)),maxshift,rotation);
    else    
        [~, tforms{m}] = xcorr2_reg(handles{1}.align,handles{m}.align,maxshift,rotation);
    end
    greenstack(:,:,m)=imtransform(handles{m}.green, tforms{m},...
        'XData',[1 size(handles{1}.green,2)],'YData',[1 size(handles{1}.green,1)]);
    redstack(:,:,m)=imtransform(handles{m}.red, tforms{m},...
        'XData',[1 size(handles{1}.red,2)],'YData',[1 size(handles{1}.red,1)]);
    ccimstack(:,:,m)=imtransform(handles{m}.ccim, tforms{m},...
        'XData',[1 size(handles{1}.ccim,2)],'YData',[1 size(handles{1}.ccim,1)])*32768; % double to uint16
end

implay(greenstack);

% 4) Create the average ccimage for ROI identification, save it to file,
% return.
ccimage=mean(ccimstack,3);

SaveWrite_MultiTiff(ccimstack,strcat(folder,'MeanCrossCorrImage.tif'));
SaveWrite_MultiTiff(greenstack,strcat(folder,'AVGgreenImage.tif'));
SaveWrite_MultiTiff(redstack,strcat(folder,'AVGredImage.tif'));
save(strcat(folder,'tforms.mat'),'tforms');

% 5) save tforms with each ccim.mat file
for m=1:length(handles)
    tform = handles{n}.tform;
    save(strcat(handles{n}.filename,'ccim.mat'), 'tform', '-append');
end

end
    
    