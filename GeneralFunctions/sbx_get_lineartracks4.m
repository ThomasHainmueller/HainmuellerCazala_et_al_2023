function [dataset] = sbx_get_lineartracks4(varargin)
% Completely revised version. Assumes the 01/17 workflow importing data
% directly from sbx files and immediately creating a 'tdata' type dataset.
%
% Extract signals from all registered image stacks using the rois provided
% from a FIJI.zip file in that folder. Categories have to be provided with
% the *.mat file corresponding to the image data
% Signaltype takes 'dGoG', 'dRoR' and 'PVcell' as arguments.
% V4 for Multiplane imaging; tforms and rois should be stored in the format
% 'commonstring'_plane1.zip etc., only the 'commonstring' is given as an
% argument!

args=struct('folder',[],'roiname','manual','tforms','tforms','minspeed',0.05,...
    'bins',.1:.05:2.1,'mode','default','channel','green','bgsubt',false); % speed in m/s
% mode: 'default','PVcells','zscored'

% Overwrite default parameters if required
for pair = reshape(varargin,2,[])
    if isfield(args,pair{1})
        args.(pair{1})=pair{2};
    else
        error('Input argument does not exist\n')
    end
end

bgsubt = args.bgsubt; %For parfor

% switch folder if necessary
if isstr(args.folder)
    folder = strrep(args.folder,'\','/');
    if ~strcmp(folder(end),'/')
        folder(end+1)='/';
    end

    cd(folder);
end

EphysFiles = ls('*.ephys');
rois = {};

if strcmp(args.channel,'green')
    chn = 1;
elseif strcmp(args.channel,'red')
    chn = 2;
else
    warning('Channel must be green or red')
    return
end

% TODO: Import all rois before parfor loop
ROIfiles = ls(['*',args.roiname,'*']);

for f = 1:size(ROIfiles,1)
    %rois{f} = ReadImageJROI(sprintf('%s\\%s_plane%d.zip',cd(), args.roiname,f));
    % Changed 171212 to allow omission of planes
    plno = str2num(ROIfiles(f,strfind(ROIfiles(f,:),'plane')+5));
    rois{plno} = ReadImageJROI(sprintf('%s\\%s_plane%d.zip',cd(), args.roiname,plno));
end

%import py.motion_corr.load_simafile

% 'Fast' Parallel loop for signal extraction
parfor n = 1:size(EphysFiles,1)
    filename = strrep(EphysFiles(n,:),'.ephys','');
    data = importfromsima(strcat(filename,'.sima')); % 5D array fr,pl,x,y,ch %EXTERNAL FKT
        
    % Step 1: Extract signals from transformed ROIs
    signals = [];
    
    % Loop over planes
    %for pl = 1:size(data,2)
    for pl = 1:length(rois)
        thistform = load(strcat(filename,'.mat'),'tform');
        
        % To maintain downward compatibility for dataset with uniform
        % tforms for all planes
        if size(thistform.tform,2) == size(data,2)
            thistform = thistform.tform{pl};
        else
            thistform = thistform.tform;
        end
        
        % MODIFIED 191013 re-implementation of dynamic background
        % subtraction
        vid = permute(squeeze(data(:,pl,:,:,chn)),[2,3,1]);
        
        if bgsubt
            vid = linewise_bg_sub2(vid);
        end
        
        if ~isempty(rois{pl})
            fprintf('%s signals extracting plane %d\n',filename,pl);
            signals = cat(2,signals, pullsignals(vid,...
                rois{pl}, thistform)); % EXTERNAL FKT % 1=signal channel for red implement later!
        else
            fprintf('%s no rois for plane %d\n',filename,pl);
        end
    end
    signalcell(n) = {signals};
    %fprintf('%s signals have been extracted\n',filename);
end

ncells = size(signalcell{1},2);

% Distribute traces from the original dataset upon the new structure
for run = 1:size(EphysFiles,1)
    filename = strrep(EphysFiles(run,:),'.ephys','');
    
    thiscat = load(filename,'category');
    thiscat = thiscat.category;
    
    % find out how many runs have been stored in this category.
    try
        catrun = length(dataset.metadata.categories{thiscat}.ft)+1;
    catch
        catrun = 1;
    end
    
    % get framerate
    info = load(filename,'info');
    info = info.info;
    
    if isempty(info.otparam)
        nplanes = 1;
    else
        nplanes = info.otparam(3);
    end
    
    framerate = info.resfreq * (2-info.scanmode) / info.config.lines / nplanes;
    
    thissignals = signalcell{run};    
    
    metadata = get_behaviourdata(strrep(EphysFiles(run,:),'.ephys',''));
    metadata = metadata(1:size(thissignals,1),:);
    
    thistform = load(strcat(filename,'.mat'),'tform');
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
            
            dataset.cells{thiscell}.image_plane = pl; %ADDED 170913
            
            %if ~args.PVcells -- replaced with args.mode 180817
            if strcmp(args.mode,'default')
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
                    SBdiscretize(thisdFoT.*thismask,...
                    dataset.metadata.categories{thiscat}.y{catrun}, args.bins);
                %dataset.cells{cell}.categories{thiscat}.dFodY{catrun} =...
                %dFoSpeed(thisdFoT,originaldata{run}.data(3,:),0:1.25E-3:4E-2);
                
            elseif strcmp(args.mode,'sheffield')
                % Sheffield et al., 2017 approach for interneurons: simply
                % divide every trace by its 8th percentile value.
                
                F0 = prctile(thissignals(:,thiscell),8);
                thisdFoT = (thissignals(:,thiscell)./F0)-1;
                sigma = nanstd(thisdFoT);
                
                dataset.cells{thiscell}.categories{thiscat}.dFoT{catrun} = thisdFoT';
                dataset.cells{thiscell}.categories{thiscat}.F0{catrun} = F0;
                dataset.cells{thiscell}.categories{thiscat}.baselineSD(catrun) = sigma;
                
                % Append ROI (only once per cell!)
                if ~isfield(dataset.cells{thiscell},'roi')
                    roi = rois{pl};
                    roi = roi{1,planecell};
                    dataset.cells{thiscell}.roi = roi;
                end
                
                % TODO: Revise this. Discretize might be inefficient!
                thismask = dataset.metadata.categories{thiscat}.moving{catrun};
                dataset.cells{thiscell}.categories{thiscat}.dFoY{catrun} =...
                    SBdiscretize(thisdFoT'.*thismask,...
                    dataset.metadata.categories{thiscat}.y{catrun}, args.bins);
                
            elseif strcmp(args.mode, 'PVcells') 
                % PVcell data extraction
                
                % PV cell approch: 1) use off-track fluorescence as
                % baseline, 2) if not possible, use F0 from last track of
                % this cat. 3) Use mean and SD of entire trace
                [thisdFoT, F0, sigma] = PVcellApproach(thissignals(:,thiscell)',...
                    dataset.metadata.categories{thiscat}.y{catrun});
                
                if isnan(F0)
                    try
                        F0 = dataset.pvcells{thiscell}.categories{thiscat}.F0{catrun-1};
                        sigma = dataset.pvcells{thiscell}.categories{thiscat}.baselineSD{catrun-1};
                    catch
                        F0 = mean(thissignals(:,thiscell));
                        sigma = std(thissignals(:,thiscell));
                    end
                    warning('cell %d category %d, run %d no baseline period available.\n',thiscell, thiscat, catrun);
                    thisdFoT = thissignals(:,thiscell)' / F0 -1;
                end
                
                dataset.pvcells{thiscell}.categories{thiscat}.dFoT{catrun} = thisdFoT;
                dataset.pvcells{thiscell}.categories{thiscat}.F0{catrun} = F0;
                dataset.pvcells{thiscell}.categories{thiscat}.baselineSD(catrun) = sigma;
                
                % Append ROI (only once per cell!)
                if ~isfield(dataset.pvcells{thiscell},'roi')
                    roi = rois{pl};
                    roi = roi{1,planecell};
                    % PARADIGM SHIFT: Struct is saved instead of mask for economic reasons
                    dataset.pvcells{thiscell}.roi = roi;
                end
                
                thismask = dataset.metadata.categories{thiscat}.moving{catrun};
                dataset.pvcells{thiscell}.categories{thiscat}.dFoY{catrun} =...
                    SBdiscretize(thisdFoT.*thismask,...
                    dataset.metadata.categories{thiscat}.y{catrun}, args.bins);          
            
            elseif strcmp(args.mode, 'zscored')
                % Z-scoring, useful for IN project.
                [thisdFoT, F0, sigma] = zscoring(thissignals(:,thiscell)');
                
                if isnan(F0)
                    try
                        F0 = dataset.cells{thiscell}.categories{thiscat}.F0{catrun-1};
                        sigma = dataset.cells{thiscell}.categories{thiscat}.baselineSD{catrun-1};
                    catch
                        F0 = mean(thissignals(:,thiscell));
                        sigma = std(thissignals(:,thiscell));
                    end
                    warning('cell %d category %d, run %d no baseline period available.\n',thiscell, thiscat, catrun);
                    thisdFoT = thissignals(:,thiscell)' / F0 -1;
                end
                
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
                
                thismask = dataset.metadata.categories{thiscat}.moving{catrun};
                dataset.cells{thiscell}.categories{thiscat}.dFoY{catrun} =...
                    SBdiscretize(thisdFoT.*thismask,...
                    dataset.metadata.categories{thiscat}.y{catrun}, args.bins);          
                
            end
        end
    end
end
end

function [thisdFoT, F0, sigma] = PVcellApproach(trace, y)
% Use off-track fluorescence as baseline
offtrack = or(y<=.05,y>=2.0);

F0 = mean(trace(offtrack));

if ~isnan(F0)
    sigma = std(trace(offtrack));
    thisdFoT = trace/F0 -1;
else
    sigma = NaN;
    thisdFoT = NaN;
end
end

function [thisdFoT, F0, sigma] = zscoring(trace)
% Use z-score trace instead of DF/F

F0 = nanmean(trace);

if ~isnan(F0)
    sigma = nanstd(trace);
    %thisdFoT = trace/F0 -1;
    thisdFoT = (trace-F0)./sigma;
else
    sigma = NaN;
    thisdFoT = NaN;
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