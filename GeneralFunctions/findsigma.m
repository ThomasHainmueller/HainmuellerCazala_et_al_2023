function [sigma, baseline] = findsigma(trace, blperiod, blfraction)
% Find the standard deviation sigma of the baseline from a calcium
% recording in vivo. Assume that the trace contains at least %blfraction of
% baseline segments of length 'blperiod' that contain no transients and
% therefore accurately represent the standard deviation of the noise.
if nargin < 3
    blperiod = 100; % 20s at 5Hz framerate.
end
if nargin < 2
    blfraction = 10;
end

for n = length(trace)-blperiod:-1:1
    allsigma(n) = nanstd(trace(n:n+blperiod));
    allbaseline(n) = nanmean(trace(n:n+blperiod));
end

sigma = prctile(allsigma,blfraction);
allsigma=sort(allsigma);
baseline = find(allsigma>=sigma,1);
baseline = allbaseline(baseline);
end