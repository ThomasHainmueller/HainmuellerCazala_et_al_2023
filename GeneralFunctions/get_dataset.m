function [masterdataset] = get_dataset(folder)
% Function to get all extracted signals from a folder and put them in a
% matrix containing frametimes, xcor, ycord, and signal of each ROI (this
% order) concatenated in one large array.

oldfolder = cd(folder);
SignalFiles = ls('*green_corr_reg_signals.csv');
masterdataset=[];

% Import data from all files in the specified folder.
for n = size(SignalFiles,1):-1:1

    filename = strrep(SignalFiles(n,:),'green_corr_reg_signals.csv','');
    SignalName = SignalFiles(n,:);
    FramesName = strcat(filename,'frametimes.csv');
    XName = strcat(filename,'xcord.csv');
    YName = strcat(filename,'ycord.csv');
    
    try
        thissignals=importdata(strcat(folder,SignalName));
        thissignals=thissignals.data(:,3:end);
        
        thisft = csvread(strcat(folder,FramesName));
        thisX = csvread(strcat(folder,XName));
        thisY = csvread(strcat(folder,YName));
        
        thisft = reshape(thisft,length(thisft),1);
        thisX = reshape(thisX,length(thisX),1);
        thisY = reshape(thisY,length(thisY),1);
        
        thisdata = cat(2,thisft,thisX,thisY,thissignals);
        masterdataset=cat(1,masterdataset,thisdata);
        
    catch
        warning([filename ' couldnt be loaded, incomplete dataset']);
    end
end


masterdataset=masterdataset';
imagesc(masterdataset(2:end,:));
k = waitforbuttonpress;

% Visual inspection of the datasets, plotX, plotY, plot signal
for n = size(masterdataset,1):-1:4
    n
    plot(masterdataset(2,:)); hold on; 
    plot(masterdataset(3,:)); hold on;
    plot(masterdataset(n,:)); hold off;
    k = waitforbuttonpress;
end


end