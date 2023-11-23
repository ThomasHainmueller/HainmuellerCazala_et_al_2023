function [SD, baseline_mean] = findsigma2(trace, eventthreshold)
% recursively determine SD for the baseline period of a calcium signal
% after removing elements that fall above the n sigma threshold.

SD = std(trace);
trace = trace-mean(trace);

while true
    baseline = trace(trace<eventthreshold*SD);
    if std(baseline) < SD
        SD = std(baseline);
    else
        baseline_mean = mean(baseline); 
        break
    end
end

end