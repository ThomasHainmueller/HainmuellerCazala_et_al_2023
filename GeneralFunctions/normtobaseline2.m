function [ntrace, F0, sigma] = normtobaseline2(trace,nsigma,iterations)
% Normalize a trace to baseline by extracting events >2 sigma from it and
% calculating the standard deviation (sigma) and baseline fluorescence F0.
if nargin < 3
    iterations = 3;
end
if nargin < 2
    nsigma = 2;
end

sigma = std(trace); % standard deviation of normalized trace
thistrace=trace;
F0 = median(thistrace);

for n = 1:iterations
    thistrace = thistrace(~(thistrace>F0+nsigma*sigma));
    sigma = std(thistrace);
    F0 = mean(thistrace);
end
ntrace = (trace-F0)/F0;
sigma=std(thistrace/F0);
end