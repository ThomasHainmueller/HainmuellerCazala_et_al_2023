function [dataset] = get_lineartracks2(folder,rois,signaltype)
% Extract signals from all registered image stacks using the rois provided
% from a FIJI.zip file in that folder and creating a 1st generation dataset
% (Should be converted to 2nd gen after categories have been set.
% Signaltype takes 'dGoG', 'dRoR' and 'PVcell' as arguments.

if nargin<3
    signaltype='dGoG';
end

cd(folder);
GreenFiles = ls('*green_corr_reg.tif');

% Number of sampling points by which xcoords are filtered for trace
% categorization
win = 25;
thisrois=ReadImageJROI([folder,rois]);

% Import data from all files in the specified folder.
for n = size(GreenFiles,1):-1:1
    filename = strrep(GreenFiles(n,:),'green_corr_reg.tif','');
    FramesName = strcat(filename,'frametimes.csv');
    XName = strcat(filename,'xcord.csv');
    YName = strcat(filename,'ycord.csv');
   
    %try
        thisft = csvread(strcat(folder,FramesName));
        thisX = csvread(strcat(folder,XName));
        thisY = csvread(strcat(folder,YName));
        % TEMPORARY BUGFIX for difference betw. existing xyt and FF
        thisft = thisft(1:length(thisX));
        
        thisft = reshape(thisft,length(thisft),1);
        thisX = reshape(thisX,length(thisX),1);
        thisY = reshape(thisY,length(thisY),1);
        thiscat = max(filter(ones(1,win)/win,1,thisX));
 
        % Categorize dataset according to x-values (maximum 2.6 for 8-bit,
        % default is ten categories.
        dataset{n}.category = numcategory(thiscat,10,2.7);
        
        % 1. Green delta F over F0 (GCaMP default)
        if strcmp(signaltype,'dGoG')
            thisgreen = ImportMultiTiff(strcat(folder,filename,'green_corr_reg.tif'));
            thisgreensignals = extract_signals(thisgreen,thisrois);
            for l=1:size(thisgreensignals,2)
                thisgreensignals(:,l)=normtobaseline2(...
                    thisgreensignals(:,l),2,3)+1.0;
            end
            thisdata = cat(2,thisft,thisX,thisY,thisgreensignals');
            dataset{n}.data=thisdata';
        end
        
        % 2. Red delta F over F0 (jRGeCO default)
        if strcmp(signaltype,'dRoR')
            thisred = ImportMultiTiff(strcat(folder,filename,'red_corr_reg.tif'));
            thisredsignals = extract_signals(thisred,thisrois);
            for l=1:size(thisredsignals,2)
                thisredsignals(:,l)=normtobaseline2(...
                    thisredsignals(:,l),2,3)+1.0;
            end
            thisdata = cat(2,thisft,thisX,thisY,thisredsignals');
            dataset{n}.data=thisdata';
        end
        
        % 3. Special mode for tdT+GCaMP, stores green dF/F plus green and
        % red raw signals
        if strcmp(signaltype,'PVcell')
            thisgreen = ImportMultiTiff(strcat(folder,filename,'green_corr_reg.tif'));
            thisgreensignals = extract_signals(thisgreen,thisrois);
            thisred = ImportMultiTiff(strcat(folder,filename,'red_corr_reg.tif'));
            thisredsignals = extract_signals(thisred,thisrois);
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
            clear thisgreensignals thisredsignals thissignals
        end       
        dataset{n}.filename=filename;
        
    %catch
        %warning([filename ' couldnt be loaded.']);
    %end
end


end