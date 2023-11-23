function [heatmap] = trial_heatmap(cellarray)
% Makes a heatmap of consecutive traces with each line representing one
% trace. Data should be organized as cell containing all the traces. Cave:
% currently this relies on equal framerates for all traces, the arrays
% returned by the autoanalysis, however, also contain the absolute
% frametimes in miliseconds.

% create a blank image with the width of the longest trace and a row for
% each trace.
rows = size(cellarray,2);
columns = max(cellfun('length',cellarray));
heatmap = nan(rows,columns);

% populate theheatmap with the traces.
for n = 1:rows
    current = cellarray{n};
    heatmap(n,1:size(current,2)) = current(1,:);
end

end