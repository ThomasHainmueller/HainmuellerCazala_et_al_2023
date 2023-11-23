% convenient function to plot all the traces in one multiple trace
% timeseries object (e.g. as returned by batchmax).
function [] = multiplot(timeser)
for c = [1:length(timeser)]
    plot(timeser(c));
    hold on;
end
end