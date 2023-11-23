function [dataset] = get_lineartracks(folder,PVcells)
% Function to get all extracted signals from a folder and put them in a
% matrix containing frametimes, xcor, ycord, and signal of each ROI (this
% order) concatenated in one large array.
if nargin<2
    PVcells=false;
end

cd(folder);
if PVcells
    SignalFiles = ls('*GCaMP_dGoR_signals.csv');
else
    SignalFiles = ls('*GCaMP_signals.csv');
end

% Number of sampling points by which xcoords are filtered for trace
% categorization
win = 25;

% Import data from all files in the specified folder.
for n = size(SignalFiles,1):-1:1
    if PVcells
        filename = strrep(SignalFiles(n,:),'GCaMP_dGoR_signals.csv','');
    else
        filename = strrep(SignalFiles(n,:),'GCaMP_signals.csv','');
    end
    SignalName = SignalFiles(n,:);
    %SignalName = strrep(SignalName,'green_corr_reg','regPV-cells');
    FramesName = strcat(filename,'frametimes.csv');
    XName = strcat(filename,'xcord.csv');
    YName = strcat(filename,'ycord.csv');
   
    %try
        thissignals=importdata(strcat(folder,SignalName));
        thissignals=thissignals.data(:,3:end);
        
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
        
        thisdata = cat(2,thisft,thisX,thisY,thissignals);
        dataset{n}.data=thisdata';
        dataset{n}.filename=filename;
        
    %catch
        %warning([filename ' couldnt be loaded, incomplete dataset']);
    %end
end

%masterdataset=masterdataset';
%imagesc(masterdataset(2:end,:));
%k = waitforbuttonpress;

% Visual inspection of the datasets, plotX, plotY, plot signal
%for n = size(masterdataset,1):-1:4
    %n
    %plot(masterdataset(2,:)); hold on; 
    %plot(masterdataset(3,:)); hold on;
    %plot(masterdataset(n,:)); hold off;
    %k = waitforbuttonpress;
%end


end