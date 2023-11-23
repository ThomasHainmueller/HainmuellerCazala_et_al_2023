function [dataset] = append_rawtracks(path,tracks,dataset)
% Append a set of raw measurements to a dataset. Conceived for
% bg-normalized Tracks from DREADD positive cells.

num_rois = floor(size(tracks,2)/4);

for n = num_rois:-1:1
    thissignals(:,n)=tracks(:,3+(n-1)*4);
end

bgval = mean(thissignals(:,end));
thissignals = thissignals(:,1:end-1)/bgval;

thisft = csvread(strcat(path,'_frametimes.csv'));
thisX = csvread(strcat(path,'_xcord.csv'));
thisY = csvread(strcat(path,'_ycord.csv'));
% TEMPORARY BUGFIX for difference betw. existing xyt and FF
thisft = thisft(1:length(thisX));

thisdata.data = cat(2,thisft',thisX',thisY',thissignals)';
thisdata.filename = path;
thisdata.category = 1;

dataset{end+1}=thisdata;
return
% Import data from all files in the specified folder.
for n = size(SignalFiles,1):-1:1

    filename = strrep(SignalFiles(n,:),'green_corr_reg_signals.csv','');
    SignalName = SignalFiles(n,:);
    %SignalName = strrep(SignalName,'green_corr_reg','regPV-cells');
    FramesName = strcat(filename,'frametimes.csv');
    XName = strcat(filename,'xcord.csv');
    YName = strcat(filename,'ycord.csv');
    
    try
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
        
    catch
        warning([filename ' couldnt be loaded, incomplete dataset']);
    end
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