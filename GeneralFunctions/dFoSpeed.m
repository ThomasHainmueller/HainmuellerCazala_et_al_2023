function dFodY = dFoSpeed(trace,y,bins)
% Calculate delta F over speed for trace and y coordinates in the range
% given by bins. Remove unlogically or <0 speeds and set to zero prior to
% calculation. Gaussian smoothed speed trace is used for calculation.
if nargin < 3
    bins = 0:1.25E-3:4E-2; % 0-32 cm/s for V2.1=4m trace and 5Hz framerate
end

dy = diff(y);
dy(end+1)=dy(end); % repeat last datapoint to get even length vector;
dy(dy<0|dy>bins(end))=0; % set 'jump' values to zero

% filter speed trace.
filter = fspecial('gaussian',[9,1],1.0);
dy = conv(dy,filter,'same');

% Generate dF over speed trace.
dFodY = SBdiscretize(trace,dy,bins);

end