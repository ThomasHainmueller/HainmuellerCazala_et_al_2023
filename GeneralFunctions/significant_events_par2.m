function [nsigma, mindur] = significant_events_par2(trace, framerate, varargin)
% Take fluorescent trace and extract the > N sigma events at a duration that
% fullfills the requirement of less than 5% negative going events.
% Procedere: Use concatenated traces from category 1 to detemine parameters
% for 0.05 chance for each cell. Then apply them for any other category.
% Don't forget to recalculate SD for each category!

args=struct('nsigma',[2,3,4],'mindur',.2,'maxdur',1,'cutoff',0.05, 'GaussFiltFreq',0,...
    'driftdur',8); % Dur: time in sec!

% Overwrite default parameters if required
for pair = reshape(varargin,2,[])
    if isfield(args,pair{1})
        args.(pair{1})=pair{2};
    else
        error('Input argument does not exist\n')
    end
end

% Gaussian filtering if required
if args.GaussFiltFreq
    gaussigma = round(framerate/args.GaussFiltFreq); % Create gauss filter for pre-processing
    gaussfilt = fspecial('gaussian', [10, 1], gaussigma);
else
    gaussfilt = 1;
end

trace = conv(trace, gaussfilt, 'same');

% Use Dombeck method for subtracting slow drift components
trace = subtract_drift(trace, round(args.driftdur*framerate), 8);

% Increase transient duration to obtain < 5% false positives
minsamples = round(args.mindur*framerate);
maxsamples = round(args.maxdur*framerate);

[SD, baseline] = findsigma2(trace, 3); % Baseline SD excluding 3sigma events
trace = trace-baseline;

for thissigma = args.nsigma
    % TODO if slow, use inner function to exclude transients of various
    % lengths.    
    for len = minsamples:maxsamples
        nPosEvents = length(find(diff(...
            transientmask2(trace, thissigma, .5, len, SD))));
        nNegEvents = length(find(diff(...
            transientmask2(-trace, thissigma, .5, len, SD))));
        if nNegEvents/nPosEvents <= args.cutoff
            nsigma = thissigma;
            mindur = len;
            return
        end
    end 
end

% In case no optimal values could be found
nsigma = args.nsigma(end);
mindur = args.maxdur;

end


function trace = subtract_drift(trace, nsamples, percentile)
% Subtract slow drift using the 8th percentile in a long time window.
% Adapted from Dombeck et al 2007.

bgtrace = zeros(length(trace),1);

for n=1:round(nsamples/2)
    bgtrace(n) = prctile(trace(1:n+round(nsamples/2)),percentile);
end

for n=round(nsamples/2)+1:length(trace)-round(nsamples/2)
    bgtrace(n) = prctile(trace(n-round(nsamples/2):n+round(nsamples/2)),...
        percentile);
end

for n=length(trace)-round(nsamples/2)+1:length(trace)
    bgtrace(n) = prctile(trace(n-round(nsamples/2):end), percentile);
end
trace = trace-bgtrace';
end