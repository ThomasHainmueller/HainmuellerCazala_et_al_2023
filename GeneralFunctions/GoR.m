% This function returns the green over red normalized matrix 
% of a given 2d line-scan matrix with the green in the first and the red channel in
% the second layer (this can be obtained from 
% the analysis window of MES by using 'matrix to variable' in the export 
% menu when both channels are selected). 
% The second return value is a vector containing the indices of
% the lines with the highest mean values in the red channel (ROI indices).

function [GoR ROI] = GoR(data, ROIsize, BGsize, F0start, F0end, tRecord)

RedMean = mean(data(:,:,2),2);
ROI = [];
BG = [];
% convert start and end times for baseline subtraction in the green channel
% from time-values (in s) to the appropriate indices of the data matrix.
F0start = round(F0start / tRecord * size(data,2));
F0end = round(F0end / tRecord * size(data,2));

% Get indices of the rows with the largest values in the red channel (ROI)
for i = 1:ROIsize
    [C I] = max(RedMean);
    ROI = [ROI I];
    RedMean(I) = NaN;
end

% Get the average of the 'BGsize' dimmest lines in the red channel as red
% Background value
for i = 1:BGsize
    [C I] = min(RedMean);
    BG = [BG C];
    RedMean(I) = NaN;
end 

BG = mean(BG);
% Background corrected red channel
red = data(:,:,2) - BG;

% Substract baselinevalue in the green channel individually for each line
for i = 1:size(data,1)
    F0 = mean(data(i,F0start:F0end,1));
    data(i,:,1) = data(i,:,1) - F0;
end

% Divide each line in the green channel by the corrosponding line in the
% red channel (corrected for the general red background BG).
for i = 1:size(data,1)
    rmean = mean(data(i,:,2))-BG;
    data(i,:,1) = data(i,:,1)/rmean;
end
GoR = data(:,:,1);
%GoR = timeseries(transpose(GoR));
%GoR.time = GoR.time*tRecord/length(GoR.time);
%image(GoR)
%GoR = data(:,:,1).*red.^-1;
end
