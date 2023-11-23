 function normalized = normalize(trace, Fbaseline)
% Convert an image stack to uint8. Min and max give the lower and upper
% boundary (in % of the total intensity distribution) of the range that is
% converted. Adjust this in order to exclude background and saturated
% pixels.

if nargin < 2
    Fbaseline = 0.05;
end

hist = sort(trace);
baseline=hist(round(Fbaseline*length(hist)));
clear hist;

normalized=trace/baseline;

end