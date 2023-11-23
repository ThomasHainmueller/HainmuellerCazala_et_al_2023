function [dataset] = get_lineartracks3(folder,rois,tforms,signaltype)
% Extract signals from all registered image stacks using the rois provided
% from a FIJI.zip file in that folder and creating a 1st generation dataset
% (Should be converted to 2nd gen after categories have been set.
% Signaltype takes 'dGoG', 'dRoR' and 'PVcell' as arguments.

if nargin<4
    signaltype='dGoG';
end

if nargin<3
    tforms = 'tforms.mat';
end

if folder
    cd(folder);
end

%GreenFiles = ls('*.align');
GreenFiles = ls('*_nonrigid.sbx');

% Number of sampling points by which xcoords are filtered for trace
% categorization
win = 25;
thisrois=ReadImageJROI([folder,rois]);

% Reload the transformation matrices which align the images.
tforms=load([folder,tforms]);
tforms=tforms.tforms;

% Import data from all files in the specified folder.
for n = size(GreenFiles,1):-1:1
    filename = strrep(GreenFiles(n,:),'_nonrigid.sbx','');
    ephysname = sprintf([filename,'.ephys']);
    
     %try
        [thisft, thisX, thisY] = get_vrposition(ephysname,3,1000); % For 1kHz standard rate.
        thiscat = max(filter(ones(1,win)/win,1,thisX));
 
        % Categorize dataset according to x-values (maximum 2.6 for 8-bit,
        % default is ten categories.
        dataset{n}.category = numcategory(thiscat,10,2.7);
        
        % Create inverse tform for this file to align rois to the image.
        thistform = fliptform(tforms{n});
        
        % 1. Green delta F over F0 (GCaMP default)
        if strcmp(signaltype,'dGoG')
            %thisgreen = ImportMultiTiff(strcat(folder,filename,'green_corr.tif'));
            %thisgreen = sbximport(filename,1,true); % Get image data;
            thisgreen = sbximport(sprintf([filename,'_nonrigid']),1,0); % For nonrigid aligned data
            [thisgreensignals, roimask] = extract_signals_transformed(...
                thisgreen,thisrois,thistform);
            dataset{n}.greenavg=mean(thisgreen,3);
            clear thisgreen
            for l=1:size(thisgreensignals,1)
                thisgreensignals(l,:)=normtobaseline2(...
                    thisgreensignals(l,:),2,3)+1.0;
            end
            thisdata = cat(2,thisft,thisX,thisY,thisgreensignals');
            dataset{n}.data=thisdata';
            clear thisgreensignals thisdata
        end
        
        % 2. Red delta F over F0 (jRGeCO default)
        if strcmp(signaltype,'dRoR')
            %thisred = ImportMultiTiff(strcat(folder,filename,'red_corr.tif'));
            thisred = sbximport(filename,2,true);
            [thisredsignals, roimask] = extract_signals_transformed(...
                thisred,thisrois,thistform);
            dataset{n}.redavg=mean(thisred,3);
            clear thisred
            for l=1:size(thisredsignals,1)
                thisredsignals(l,:)=normtobaseline2(...
                    thisredsignals(l,:),2,3)+1.0;
            end
            thisdata = cat(2,thisft,thisX,thisY,thisredsignals');
            dataset{n}.data=thisdata';
            clear thisredsignals thisdata
        end
        
        % 3. Special mode for tdT+GCaMP, stores green dF/F plus green and
        % red raw signals
        if strcmp(signaltype,'PVcell')
            thisgreen = sbximport(filename,1,true);
            [thisgreensignals, roimask] = extract_signals_transformed(...
                thisgreen,thisrois,thistform);            
            thisgreen = sbximport(filename,2,true);
            thisredsignals = extract_signals_transformed(...
                thisred,thisrois,thistform);
            dataset{n}.greenavg=mean(thisgreen,3);
            dataset{n}.redavg=mean(thisred,3);
            clear thisgreen thisred
            for l=size(thisgreensignals,1):-1:1
                % Conventional baseline to 10% of lowest intensity values
                thissignals(l,:)=normtobaseline(thisgreensignals(l,:),0.1)+1.0;
            end
            thisdata = cat(2,thisft,thisX,thisY,thissignals');
            thisgreenraw = cat(2,thisft,thisX,thisY,thisgreensignals');
            thisredraw = cat(2,thisft,thisX,thisY,thisredsignals');
            dataset{n}.data=thisdata';
            dataset{n}.greenraw=thisgreenraw';
            dataset{n}.redraw=thisredraw';
            clear thisgreensignals thisredsignals thissignals thisdata
        end       
        
        % 4. Special mode for tdT+GCaMP, stores green dF/F plus green and
        % red raw signals for independent ROI in each run.
        if strcmp(signaltype,'PVcell2')
            if mod(length(thisrois),length(GreenFiles))~=0
                warning('Number of ROIs must be a multiple of the number of runs!');
                return
            end
            ncells=length(thisrois)/length(GreenFiles);
            thisgreen = sbximport(filename,1,true);
            % Get appropriate ROIs for the present run
            %%%thisrunsrois=thisrois([1:ncells]*n);
            thisrunsrois=thisrois((([1:ncells]-1)*length(GreenFiles))+n);
            [thisgreensignals, roimask] = extract_signals_transformed(...
                thisgreen,thisrunsrois,thistform);            
            thisgreen = sbximport(filename,2,true);
            thisredsignals = extract_signals_transformed(...
                thisred,thisrois,thistform);
            dataset{n}.greenavg=mean(thisgreen,3);
            dataset{n}.redavg=mean(thisred,3);
            clear thisgreen thisred
            for l=size(thisgreensignals,1):-1:1
                % Conventional baseline to 10% of lowest intensity values
                thissignals(l,:)=normtobaseline(thisgreensignals(l,:),0.1)+1.0;
            end
            thisdata = cat(2,thisft,thisX,thisY,thissignals');
            thisgreenraw = cat(2,thisft,thisX,thisY,thisgreensignals');
            thisredraw = cat(2,thisft,thisX,thisY,thisredsignals');
            dataset{n}.data=thisdata';
            dataset{n}.greenraw=thisgreenraw';
            dataset{n}.redraw=thisredraw';
            clear thisgreensignals thisredsignals thissignals thisdata
        end       
        
        dataset{n}.filename=filename;
        
        % Append ROI overview to file.
        dataset{n}.roimask=max(roimask,[],3);  
%     catch
%         warning([filename ' couldnt be loaded.']);
%     end
warning([filename ' has been processed.']);
end
end