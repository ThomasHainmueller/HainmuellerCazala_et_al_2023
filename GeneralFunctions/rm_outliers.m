function trace = rm_outliers(trace, Fout, replace, direction)
% Remove the Fout % of datapoints from a trace and set to 'replace' value;
% direction can be 'positive', 'negative' or 'both' to indicate if plus or
% minus signed outliers should be considered. Usefull for removing
% artificial 'reset' signals from speed traces.

if nargin < 4
    direction = 'both';
end

if nargin < 3
    replace = 0;
end

if nargin < 2
    Fout = 0.05;
end

hist = sort(trace);
lowest = hist(round(Fout*length(trace)));
highest = hist(round((1-Fout)*length(trace)));

if strcmp(direction,'positive')
    trace(trace>highest) = replace;
elseif strcmp(direction,'negative')
    trace(trace<lowest) = replace;
elseif strcmp(direction,'both')
    trace(trace>highest) = replace;
    trace(trace<lowest) = replace;
end
end