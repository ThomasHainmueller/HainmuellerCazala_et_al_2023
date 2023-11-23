function corrected = remove_fluctuations(trace, window, percentile)
% Takes a trace of continuous data and subtracts from each datapoint the
% given percentile of the values in a given window around the datapoint.
% Thereby removes slow (artificial) fluctuations from traces containing
% fast signals.
corrected = [];
w = floor(window/2);
last = length(trace);

% treating ends
for n = last:-1:last-w
    corrected(n) = trace(n)-prctile(trace(last-w:last),percentile);
end

for n = 1:w
    corrected(n) = trace(n)-prctile(trace(1:w),percentile);
end

for n = w+1:last-w
    corrected(n) = trace(n)-prctile(trace(n-w:n+w),percentile);
end
end