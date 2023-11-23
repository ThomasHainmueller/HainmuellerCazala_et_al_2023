function transients = transientmask2(trace, nSigmaStart, nSigmaEnd, mindur, SD)
% Create mask for transients exceeding the nSigmaStart level lasting until
% crossing nSigmaEnd again with a specified minimum duration. Adapted from
% Dombeck et al., 2007. Mindur is in SAMPLING POINTS!
if nargin<5
    SD = findsigma2(trace, nSigmaStart); % find baseline sigma
end

% Find transient start or end of trace, whatever is next
if find(trace >= nSigmaStart*SD, 1, 'first')
    first = find(trace >= nSigmaStart*SD, 1, 'first');
else
    transients = false(1,length(trace));
    return
end

% Find end of transient or end of trace, whatever next
if find(trace(first:end) <= nSigmaEnd*SD, 1, 'first')
    dur = find(trace(1,first:end) <= nSigmaEnd*SD, 1, 'first');
else
    dur = length(trace);
end

% RECURSION PART:
% Check for mindur criterion, set zero if testing is not requested
if dur >= mindur
    % Case: transient is longer than mindur (or any length if mindur=0)
    transients = [false(1,first), true(1,dur),...
        transientmask2(trace(first+dur+1:end), nSigmaStart, nSigmaEnd, mindur, SD)];  
else
    % Case: transient is shorter than mindur - do not include
    transients = [false(1,first+dur),...
        transientmask2(trace(first+dur+1:end), nSigmaStart, nSigmaEnd, mindur, SD)];
end
end