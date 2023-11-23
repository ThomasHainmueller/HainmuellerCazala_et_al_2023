function [dataset] = tiff_get_lineartracks4(folder,rois,tforms,nplanes,signaltype,separatefolders)
% Extract signals from all registered image stacks using the rois provided
% from a FIJI.zip file in that folder and creating a 1st generation dataset
% (Should be converted to 2nd gen after categories have been set.
% Signaltype takes 'dGoG', 'dRoR' and 'PVcell' as arguments.
% V4 for Multiplane imaging; tforms and rois should be stored in the format
% 'commonstring'_plane1.zip etc., only the 'commonstring' is given as an
% argument!
if nargin<6
    separatefolders=true; % Files for individual planes stay in separate folders
end
if nargin<5
    signaltype='dGoG';
end

% format folder string as necessarry
folder = strrep(folder,'\','/');
if ~strcmp(folder(end),'/')
    folder(end+1)='/';
end

cd(folder);
EphysFiles = ls('*.ephys');

% Number of sampling points to filter xcoords for categorization
win = 25;

% Loop through recording filesets in folder.
for n = size(EphysFiles,1):-1:1
    ephysname = EphysFiles(n,:);
    filename = strrep(EphysFiles(n,:),'.ephys','');
    
    try
        [thisft, thisX, thisY] = get_vrposition(ephysname,3,1000,'tiff','bidirectional'); % For 1kHz standard rate.
        nframes = length(thisft); % Number of total frames. Note: Each plane has only nframes/nplanes images!
        
        % TEMPORARY BUGFIX for difference betw. existing xyt and FF
        thiscat = max(filter(ones(1,win)/win,1,thisX));
        
        % Categorize dataset according to x-values (maximum 2.6 for 8-bit,
        % default is ten categories.
        dataset{n}.category = numcategory(thiscat,10,2.7);
        dataset{n}.data = cat(1,thisft',thisX',thisY');
        dataset{n}.filename=filename;
    
        for pl=1:nplanes
            if separatefolders
                thisrois=ReadImageJROI(sprintf([folder,'/plane%d/',rois,'_plane%d.zip'],pl,pl));        
                % Reload the transformation matrices which align the images.
                thistforms=load(sprintf([folder,'/plane%d/',tforms,'_plane%d.mat'],pl,pl));                
            else
                thisrois=ReadImageJROI(sprintf([folder,rois,'_plane%d.zip'],pl));
                thistforms=load(sprintf([folder,tforms,'_plane%d.mat'],pl));  
            end
            
            thistforms=thistforms.tforms;

            % Create inverse tform for this file to align rois to the image.
            thistform = fliptform(thistforms{n});

            % 1. Green delta F over F0 (GCaMP default)
            if strcmp(signaltype,'dGoG')
                %thisgreen = ImportMultiTiff(strcat(folder,filename,'_green_corr.tif'));
                if separatefolders
                    thisgreen = ImportMultiTiff(sprintf([folder,'/plane%d/',filename,'_plane%d_green_corr.tif'],pl,pl));
                else
                    thisgreen = ImportMultiTiff(sprintf([folder,filename,'_plane%d_green_corr.tif'],pl));
                end
                [thisgreensignals, roimask] = extract_signals_transformed(...
                    thisgreen,thisrois,thistform);
                dataset{n}.greenavg{pl}=mean(thisgreen,3);
                clear thisgreen
                for l=1:size(thisgreensignals,1)
                    thisgreensignals(l,:)=normtobaseline2(...
                        thisgreensignals(l,:),2,3)+1.0;
                end
                %thisdata = cat(2,thisft,thisX,thisY,thisgreensignals');
                thisgreensignals=interp1(pl:nplanes:nframes,thisgreensignals',1:nframes,'linear',1.0); % Interpolate missing values to match the three planes to their respective x/y/ft data
                dataset{n}.data=cat(1,dataset{n}.data,thisgreensignals');
                clear thisgreensignals thisdata
            end

%             % 2. Red delta F over F0 (jRGeCO default)
%             if strcmp(signaltype,'dRoR')
%                 thisred = ImportMultiTiff(strcat(folder,filename,'_red_corr.tif'));
%                 [thisredsignals, roimask] = extract_signals_transformed(...
%                     thisred,thisrois,thistform);
%                 dataset{n}.redavg=mean(thisred,3);
%                 clear thisred
%                 for l=1:size(thisredsignals,1)
%                     thisredsignals(l,:)=normtobaseline2(...
%                         thisredsignals(l,:),2,3)+1.0;
%                 end
%                 thisdata = cat(2,thisft,thisX,thisY,thisredsignals');
%                 dataset{n}.data=thisdata';
%                 clear thisredsignals thisdata
%             end
% 
%             % 3. Special mode for tdT+GCaMP, stores green dF/F plus green and
%             % red raw signals
%             if strcmp(signaltype,'PVcell')
%                 thisgreen = ImportMultiTiff(strcat(folder,filename,'_green_corr.tif'));
%                 [thisgreensignals, roimask] = extract_signals_transformed(...
%                     thisgreen,thisrois,thistform);            
%                 thisred = ImportMultiTiff(strcat(folder,filename,'_red_corr.tif'));
%                 thisredsignals = extract_signals_transformed(...
%                     thisred,thisrois,thistform);
%                 dataset{n}.greenavg=mean(thisgreen,3);
%                 dataset{n}.redavg=mean(thisred,3);
%                 clear thisgreen thisred
%                 for l=size(thisgreensignals,1):-1:1
%                     % Conventional baseline to 10% of lowest intensity values
%                     thissignals(l,:)=normtobaseline(thisgreensignals(l,:),0.1)+1.0;
%                 end
%                 thisdata = cat(2,thisft,thisX,thisY,thissignals');
%                 thisgreenraw = cat(2,thisft,thisX,thisY,thisgreensignals');
%                 thisredraw = cat(2,thisft,thisX,thisY,thisredsignals');
%                 dataset{n}.data=thisdata';
%                 dataset{n}.greenraw=thisgreenraw';
%                 dataset{n}.redraw=thisredraw';
%                 clear thisgreensignals thisredsignals thissignals thisdata
%             end       
% 
%             % 4. Special mode for tdT+GCaMP, stores green dF/F plus green and
%             % red raw signals for independent ROI in each run.
%             if strcmp(signaltype,'PVcell2')
%                 if mod(length(thisrois),length(EphysFiles))~=0
%                     warning('Number of ROIs must be a multiple of the number of runs!');
%                     return
%                 end
%                 ncells=length(thisrois)/length(EphysFiles);
%                 thisgreen = ImportMultiTiff(strcat(folder,filename,'_green_corr.tif'));
%                 % Get appropriate ROIs for the present run
%                 %%%thisrunsrois=thisrois([1:ncells]*n);
%                 thisrunsrois=thisrois((([1:ncells]-1)*length(EphysFiles))+n);
%                 [thisgreensignals, roimask] = extract_signals_transformed(...
%                     thisgreen,thisrunsrois,thistform);            
%                 thisred = ImportMultiTiff(strcat(folder,filename,'_red_corr.tif'));
%                 thisredsignals = extract_signals_transformed(...
%                     thisred,thisrois,thistform);
%                 dataset{n}.greenavg=mean(thisgreen,3);
%                 dataset{n}.redavg=mean(thisred,3);
%                 clear thisgreen thisred
%                 for l=size(thisgreensignals,1):-1:1
%                     % Conventional baseline to 10% of lowest intensity values
%                     thissignals(l,:)=normtobaseline(thisgreensignals(l,:),0.1)+1.0;
%                 end
%                 thisdata = cat(2,thisft,thisX,thisY,thissignals');
%                 thisgreenraw = cat(2,thisft,thisX,thisY,thisgreensignals');
%                 thisredraw = cat(2,thisft,thisX,thisY,thisredsignals');
%                 dataset{n}.data=thisdata';
%                 dataset{n}.greenraw=thisgreenraw';
%                 dataset{n}.redraw=thisredraw';
%                 clear thisgreensignals thisredsignals thissignals thisdata
%             end       

            % Append ROI overview to file.
            dataset{n}.roimask{pl}=max(roimask,[],3);
            fprintf([filename,'_plane%d has been processed.\n'],pl);
        end
    catch
        warning([filename ' couldnt be loaded.']);
    end
end

end
