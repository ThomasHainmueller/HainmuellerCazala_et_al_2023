function [greenavg,redavg] = register_images(folder, maxshift, RedCorr)
% Analysis routine that accepts G/R datasets, aligns them according to the
% red channel. Then computes the Cross-correlation image of the green image
% and finally writes the registered green, red and cross corrs to .tif. In
% addition, it generates and saves a mean cross-corr image of all
% registered datasets which can be used for ROI selection in FIJI.

if nargin < 4
    % Use the mean of all images registered to this point instead of the
    % first image for registration.
    usemean = false;
end

if nargin < 3
    RedCorr = true;
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
    RedName = strcat(filename,'red_corr.tif');
    
    try
        % Load background channel. Named 'red' for historical reasons even
        % if it is in fact the mean of the green channel!
        if RedCorr
            thishandle.red = mean(ImportMultiTiff(strcat(folder,RedName)),3);
            thishandle.red=BGsubt(thishandle.red,0.001);
        elseif ~RedCorr
            thishandle.red = mean(ImportMultiTiff(strcat(folder,GreenName)),3);
            %Rolling Ball BG-subtraction:
            bg = imopen(thishandle.red,strel('disk',2));
            thishandle.red = thishandle.red-bg;
        end
        thishandle.gname=GreenName;
        thishandle.rname=RedName;
        handles{end+1}=thishandle;
    catch
        warning([filename ' couldnt be loaded, incomplete dataset']);
    end
end

% 2) Register all image stacks; Make a stack of red reference images after
% registration (for finding the image area common to all datasets and
% append the transformation information for registration to the handle of
% all datasets.

% Create registration routine
%[optimizer, metric]=imregconfig('monomodal');
%optimizer.MaximumIterations = 1e5;
%optimizer.MaximumStepLength = 0.8;
redrefstack=handles{1}.red;


% register each red reference image
for m=2:length(handles)
    reference = handles{1}.red;
    
    % Make averaged reference image if required.
    %if usemean
        %for im = 2:m-1
            %reference = cat(3,reference,handles{im}.red);
        %end
        %reference = mean(reference,3);
    %end
    
    [thisredreg, thistform] = xcorr2_reg(reference,handles{m}.red,maxshift);
    handles{m}.tform = thistform;
    
    redrefstack=cat(3,redrefstack,thisredreg);
    clear thisredreg thistform;
end

implay(redrefstack);
[xbound,ybound]=imcrop_coords(redrefstack);

% 3) Register and crop each dataset to the common area for all datasets,
% write the modified data to tiff. Generate ccimages for each dataset,
% write them to tiff and append them to the ccimstack.
greenavg = [];
redavg = [];

for m=length(handles):-1:1
    thisgreen = ImportMultiTiff(strcat(folder,handles{m}.gname));
    thisred = ImportMultiTiff(strcat(folder,handles{m}.rname));
    
    % For all Datasets except the first: register green and red corr
    % datasets and crop them to the common boundaries.
    if m ~= 1
        fixed = handles{1}.red;
        thisred = imtransform(thisred,handles{m}.tform,'XData',[1 size(fixed,2)],'YData',[1 size(fixed,1)]);
        thisgreen = imtransform(thisgreen,handles{m}.tform,'XData',[1 size(fixed,2)],'YData',[1 size(fixed,1)]);
    end
    
    thisred = thisred(xbound,ybound,:);
    thisgreen=thisgreen(xbound,ybound,:);
    
    % Optional background subtraction. Usefull if images from datasets are
    % to be concatenated to masterstacks later on.
    %thisgreen=BGsubt(thisgreen,0.001);
    %thisred=BGsubt(thisred,0.001);
    
    greenavg = cat(3,greenavg,mean(thisgreen,3));
    redavg = cat(3,redavg,mean(thisred,3));

    % write tif files
    gregname=strrep(handles{m}.gname,'corr.tif','corr_reg.tif');
    rregname=strrep(handles{m}.rname,'corr.tif','corr_reg.tif');

    SaveWrite_MultiTiff(thisgreen,strcat(folder,gregname));
    SaveWrite_MultiTiff(thisred,strcat(folder,rregname));
    
    clear thisgreen thisred thisccim;
end

% 4) Create the average ccimage for ROI identification, save it to file,
% return.

SaveWrite_MultiTiff(uint16(greenavg),strcat(folder,'AVGgreenStack.tif'));
SaveWrite_MultiTiff(uint16(redavg),strcat(folder,'AVGredStack.tif'));

end
    
    