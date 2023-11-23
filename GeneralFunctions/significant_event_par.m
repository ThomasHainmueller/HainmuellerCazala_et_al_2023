function [nsigma, dT, pval] = significant_event_par(trace,sigmarange,timerange, cutoff)
% Determine the minimum values of sigma and dT from a given parameter range
% that yield a detection rate of less than 'cutoff' false positives and
% at the maximum rate of detection for positive events.

if nargin < 4
    cutoff = 0.05;
end

% Find the standard deviation of the baseline for all categories
blperiod = 100; % 20s at 5Hz framerate.
blfraction = 10;
sigma = findsigma(trace,blperiod,blfraction);

for s = 1:length(sigmarange)
    nsigma = sigmarange(s);
    for t = 1:length(timerange)
        dT = timerange(t);
        % Loop over trace, set positive or negative transients as positive
        % every time the trace undercuts or exceeds nsigma for dT.
        [positive,negative] = deal(false(1,length(trace)));
        n = 1;
        while n<=length(trace)-dT
            if trace(n:n+dT)>nsigma*sigma
                positive(n)=1;
                n=n+dT+1;
            elseif trace(n:n+dT)<-nsigma*sigma
                negative(n)=1;
                n=n+dT+1;
            else
                n=n+1;
            end
        end
        error_rate(s,t)=length(find(negative))/length(find(positive));
        events_detected(s,t)=length(find(positive));
    end
end
correct_mask=error_rate<cutoff;
if ~any(correct_mask>0)
    [~,index] = min(error_rate(:));
    [bestsigma,bestdT] = ind2sub(size(error_rate),index);
else
    events_detected = events_detected.*correct_mask;
    [~,index] = max(events_detected(:));
    [bestsigma,bestdT] = ind2sub(size(events_detected),index);
end
nsigma=sigmarange(bestsigma);
dT = timerange(bestdT);
pval = error_rate(bestsigma,bestdT);
end