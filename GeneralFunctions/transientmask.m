function mask = transientmask(trace, dTmin, nsigmarise, nsigmadecay, sigma)
% Make a mask that identifies the significant transients in a calcium
% recoridng that cross the nsigmarise and last until crossing nsigmadecay.
% Adapted from Dombeck et al., 2007.
if nargin < 5
    sigma = findsigma(trace,100,0.1);
end
if nargin < 4
    nsigmadecay = 0.5;
end
if nargin < 3
    nsigmarise = 2.0;
end
if nargin < 2
    dTmin = 0;
end

mask = false(1,length(trace));

for n = 1:length(trace)-dTmin
    if trace(n) > nsigmarise*sigma
        if trace(n+dTmin)>nsigmarise*sigma
            last = find(trace(n:end)<nsigmadecay*sigma,1);
            if ~isempty(last) && n+last<length(trace)
                trace(n:n+last)=0;
                mask(n:n+last)=1;
            else
                mask(n:length(trace))=1;
            end
        end
    end
end
end