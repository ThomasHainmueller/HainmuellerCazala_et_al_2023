function [dataset] = sbx_get_lineartracks4(varargin)
% Completely revised version. Assumes the 01/17 workflow importing data
% directly from tiff files and immediately creating a 'tdata' type dataset.
% This is conceived for the remaining femtonics data.
%
% Extract signals from all registered image stacks using the rois provided
% from a FIJI.zip file in that folder. Categories have to be provided with
% the *ccim.mat file corresponding to the image data
% Signaltype takes 'dGoG', 'dRoR' and 'PVcell' as arguments.
% V4 for Multiplane imaging; tforms and rois should be stored in the format
% 'commonstring'_plane1.zip etc., only the 'commonstring' is given as an
% argument!

args=struct('folder',[],'roiname','manual','tforms','tforms','minspeed',0.05,...
    'bins',.1:.05:2.1); % speed in m/s

% Overwrite default parameters if required
for pair = reshape(varargin,2,[])
    if isfield(args,pair{1})
        args.(pair{1})=pair{2};
    else
        error('Input argument does not exist\n')
    end
end

% switch folder if necessary
if isstr(args.folder)
    folder = strrep(args.folder,'\','/');
    if ~strcmp(folder(end),'/')
        folder(end+1)='/';
    end

    cd(folder);
end

TiffFiles = ls('*green_corr.tif');
rois = {};

% TODO: Import all rois before parfor loop
ROIfiles = ls(['*',args.roiname,'*']);
for f = 1:size(ROIfiles,1)
    rois{f} = ReadImageJROI(sprintf('%s\\%s_plane%d.zip',cd(),args.roiname,f));
end

%import py.motion_corr.load_simafile

% 'Fast' Parallel loop for signal extraction
parfor n = 1:size(TiffFiles,1)
    filename = strrep(TiffFiles(n,:),'_green_corr.tif','');
    data = [];
    
    % Get data in format usually broadcasted by sbx files
    data(:,1,:,:,1) = permute(ImportMultiTiff(strcat(filename,'_green_corr.tif')),[3,1,2]);
    data(:,1,:,:,2) = permute(ImportMultiTiff(strcat(filename,'_red_corr.tif')),[3,1,2]);
        
    % Step 1: Extract signals from transformed ROIs
    signals = [];
    
    % Loop over planes
    for pl = 1:size(data,2)
        thistform = load(strcat(filename,'_ccim.mat'),'tform');
        thistform = thistform.tform;
        
        signals = cat(2,signals, pullsignals(permute(squeeze(data(:,pl,:,:,1)),[2,3,1]),...
            rois{pl}, thistform)); % EXTERNAL FKT % 1=signal channel for red implement later!
    end
    signalcell(n) = {signals};
    fprintf('%s signals have been extracted\n',filename);
end

ncells = size(signalcell{1},2);

% Distribute traces from the original dataset upon the new structure
for run = 1:size(TiffFiles,1)
    filename = strrep(TiffFiles(run,:),'_green_corr.tif','');
    
    thiscat = load(strcat(filename,'_ccim.mat'),'category');
    thiscat = thiscat.category;
    
    % find out how many runs have been stored in this category.
    try
        catrun = length(dataset.metadata.categories{thiscat}.ft)+1;
    catch
        catrun = 1;
    end
    
    % get framerate
    framerate = 1000/mean(diff(csvread(strcat(...
        filename,'_frametimes.csv')))); % In Hz
    
    thissignals = signalcell{run};

    metadata(:,1) = csvread(strcat(filename,'_xcord.csv'))';
    metadata(:,2) = csvread(strcat(filename,'_ycord.csv'))';
    % Fill lick- and pupil-data with NaNs;
    [metadata(:,3),metadata(:,4)] = deal(nan(size(metadata(:,1),1),1));
    
    thistform = load(strcat(filename,'_ccim.mat'),'tform');
    thistform = thistform.tform;        
    
    dataset.metadata.categories{thiscat}.acquisition_rate(catrun) = framerate;
    dataset.metadata.categories{thiscat}.ft{catrun} = ...
        [1:size(metadata,1)]*1000/framerate; % in milliseconds
    dataset.metadata.categories{thiscat}.x{catrun} = metadata(:,1);
    dataset.metadata.categories{thiscat}.y{catrun} = metadata(:,2);
    dataset.metadata.categories{thiscat}.licks{catrun} = metadata(:,3);
    dataset.metadata.categories{thiscat}.pupil_area{catrun} = metadata(:,4);
    dataset.metadata.categories{thiscat}.filename{catrun} = filename;
    dataset.metadata.categories{thiscat}.moving{catrun}=...
         movingmask2(metadata(:,2), framerate, args.minspeed); %Removed Matrix inversion 170519
    dataset.metadata.categories{thiscat}.tforms{catrun} = thistform;
    
    % DEPRECATED: Background images (ca. 2-3 GB total) are stored with
    % individual .mat files now as .images porperty.
    
    thiscell = 0;   
    for pl = 1:length(rois)    
        for planecell = 1:length(rois{pl})
            thiscell = thiscell +1;
            
            % TODO: Reimplement PV cells, impl. new baseline determination!
            [thisdFoT,F0,sigma] = normtobaseline2(thissignals(:,thiscell)',2,5); %170528 Revised vector orientation to 1xn
            dataset.cells{thiscell}.categories{thiscat}.dFoT{catrun} = thisdFoT;
            dataset.cells{thiscell}.categories{thiscat}.F0{catrun} = F0;
            dataset.cells{thiscell}.categories{thiscat}.baselineSD(catrun) = sigma;
            
            % Append ROI (only once per cell!)
            if ~isfield(dataset.cells{thiscell},'roi')
                roi = rois{pl};
                roi = roi{1,planecell};
                % PARADIGM SHIFT: Struct is saved instead of mask for economic reasons
                dataset.cells{thiscell}.roi = roi;
            end
            
            % TODO: Revise this. Discretize might be inefficient!
            thismask = dataset.metadata.categories{thiscat}.moving{catrun};
            dataset.cells{thiscell}.categories{thiscat}.dFoY{catrun} =...
                discretize(thisdFoT.*thismask,...
                dataset.metadata.categories{thiscat}.y{catrun}, args.bins);
            %dataset.cells{cell}.categories{thiscat}.dFodY{catrun} =...
                %dFoSpeed(thisdFoT,originaldata{run}.data(3,:),0:1.25E-3:4E-2);         
        end
    end
    clear metadata
end
end

%     %for cell = ncells:-1:1
%     % PARADIGM SHIFT: Struct is saved instead of mask for economic reasons
%     %dataset.metadata.categories{thiscat}.roi{catrun} = rois{run};
%     %dataset.metadata.categories{thiscat}.tform{catrun} = tforms{run};
%         % Calculate baseline and bl sigma of the trace without transients
%         if PVcells
%             thisdFoT=originaldata{run}.data(cell+3,:)-1.0;
%             thisGraw=originaldata{run}.greenraw(cell+3,:);
%             thisRraw=originaldata{run}.redraw(cell+3,:);
%             dataset.cells{cell}.categories{thiscat}.greenraw{catrun}=thisGraw;
%             dataset.cells{cell}.categories{thiscat}.redraw{catrun}=thisRraw;
%             dataset.cells{cell}.categories{thiscat}.meanGoR(catrun)=...
%                 mean(thisGraw)/mean(thisRraw);
%             %dataset.cells{cell}.categories{thiscat}.dFodY{catrun} =...
%                 %dFoSpeed(thisdFoT,originaldata{run}.data(3,:),0:1.25E-3:4E-2);
%         else  
%             [thisdFoT,F0,sigma] = normtobaseline2(originaldata{run}.data(cell+3,:),2,5);
%             dataset.cells{cell}.categories{thiscat}.baseline(catrun) = F0;
%             dataset.cells{cell}.categories{thiscat}.baselineSD(catrun) = sigma;
%         end
%         thismask = dataset.metadata.categories{thiscat}.moving{catrun};
%         dataset.cells{cell}.categories{thiscat}.dFoT{catrun} = thisdFoT;
%         dataset.cells{cell}.categories{thiscat}.dFoY{catrun} =...
%             discretize(thisdFoT.*thismask,originaldata{run}.data(3,:),bins);
%         % Calculate dFodY (running speed), default 0-32cm/s range is
%         % assumed; For the entire (non-trace based) Calcium signal.
%         dataset.cells{cell}.categories{thiscat}.dFodY{catrun} =...
%             dFoSpeed(thisdFoT,originaldata{run}.data(3,:),0:1.25E-3:4E-2);
%     end
% end
% Loop through recording filesets in folder.
% for n = size(EphysFiles,1):-1:1
%     ephysname = EphysFiles(n,:);
%     filename = strrep(EphysFiles(n,:),'.ephys','');
%     
%     try
%         [thisft, thisX, thisY] = get_vrposition(ephysname,3,1000,'tiff','bidirectional'); % For 1kHz standard rate.
%         nframes = length(thisft); % Number of total frames. Note: Each plane has only nframes/nplanes images!
%         
%         % TEMPORARY BUGFIX for difference betw. existing xyt and FF
%         thiscat = max(filter(ones(1,win)/win,1,thisX));
%         
%         % Categorize dataset according to x-values (maximum 2.6 for 8-bit,
%         % default is ten categories.
%         dataset{n}.category = numcategory(thiscat,10,2.7);
%         dataset{n}.data = cat(1,thisft',thisX',thisY');
%         dataset{n}.filename=filename;
%     
%         for pl=1:nplanes
%             if separatefolders
%                 thisrois=ReadImageJROI(sprintf([folder,'/plane%d/',rois,'_plane%d.zip'],pl,pl));        
%                 % Reload the transformation matrices which align the images.
%                 thistforms=load(sprintf([folder,'/plane%d/',tforms,'_plane%d.mat'],pl,pl));                
%             else
%                 thisrois=ReadImageJROI(sprintf([folder,rois,'_plane%d.zip'],pl));
%                 thistforms=load(sprintf([folder,tforms,'_plane%d.mat'],pl));  
%             end
%             
%             thistforms=thistforms.tforms;
% 
%             % Create inverse tform for this file to align rois to the image.
%             thistform = fliptform(thistforms{n});
% 
%             % 1. Green delta F over F0 (GCaMP default)
%             if strcmp(signaltype,'dGoG')
%                 %thisgreen = ImportMultiTiff(strcat(folder,filename,'_green_corr.tif'));
%                 if separatefolders
%                     thisgreen = ImportMultiTiff(sprintf([folder,'/plane%d/',filename,'_plane%d_green_corr.tif'],pl,pl));
%                 else
%                     thisgreen = ImportMultiTiff(sprintf([folder,filename,'_plane%d_green_corr.tif'],pl));
%                 end
%                 [thisgreensignals, roimask] = extract_signals_transformed(...
%                     thisgreen,thisrois,thistform);
%                 dataset{n}.greenavg{pl}=mean(thisgreen,3);
%                 clear thisgreen
%                 for l=1:size(thisgreensignals,1)
%                     thisgreensignals(l,:)=normtobaseline2(...
%                         thisgreensignals(l,:),2,3)+1.0;
%                 end
%                 %thisdata = cat(2,thisft,thisX,thisY,thisgreensignals');
%                 thisgreensignals=interp1(pl:nplanes:nframes,thisgreensignals',1:nframes,'linear',1.0); % Interpolate missing values to match the three planes to their respective x/y/ft data
%                 dataset{n}.data=cat(1,dataset{n}.data,thisgreensignals');
%                 clear thisgreensignals thisdata
%             end   
% 
%             % Append ROI overview to file.
%             dataset{n}.roimask{pl}=max(roimask,[],3);
%             fprintf([filename,'_plane%d has been processed.\n'],pl);
%         end
%     catch
%         warning([filename ' couldnt be loaded.']);
%     end
% end