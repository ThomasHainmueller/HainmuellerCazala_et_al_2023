function data = significant_events2(data, varargin)
% Master function to handle dataset, distribute recording traces and store
% the generated masks for transients.

args=struct('nsigma',2:5,'mindur',.1,'maxdur',1,'cutoff',0.05, 'GaussFiltFreq',1,...
    'driftdur',8,'framerate','default'); % Dur: time in sec!
%args=struct('nsigma',3,'mindur',.2,'maxdur',5,'cutoff',0.05, 'GaussFiltFreq',1,...
    %'driftdur',8,'framerate','default'); % Dur: time in sec!-Adapted 181021

% Overwrite default parameters if required
for pair = reshape(varargin,2,[])
    if isfield(args,pair{1})
        args.(pair{1})=pair{2};
    else
        error('Input argument does not exist\n')
    end
end

% Can be specified for older datasets w.o. stored framerates
if strcmp(args.framerate,'default')
    args.framerate = data.metadata.categories{1}.acquisition_rate(1);
end

% Loop over cells in dataset
for n = 1:length(data.cells)
    % Collect and collate all calcium transients from this cell
    thisdFoT = [];
    for c = 1:length(data.cells{n}.categories)
        traces = deal(data.cells{n}.categories{c}.dFoT);
        traces = cat(2,traces{:});
        traces = traces-median(traces);
        thisdFoT = cat(2,thisdFoT,traces);
        clear traces
    end
    
    
    % Find the values for nsigma and dT that yield the highest detection
    % rate within the cutoff tolerance of error rate.
    [nsigma, mindur]=significant_events_par2(thisdFoT, args.framerate,...
        'nsigma',args.nsigma,'mindur',args.mindur,'maxdur',args.maxdur,'cutoff',args.cutoff,...
        'GaussFiltFreq',args.GaussFiltFreq,'driftdur',15);
    data.cells{n}.transient_sigma = nsigma;
    data.cells{n}.transient_mindur = mindur;
    
    % Create masks for all traces in all categories and store them with the
    % dataset.
    for c = 1:length(data.cells{n}.categories)
        %data.cells{n}.categories{c}.baselinesigma = sigma(c);
        for tr = 1:length(data.cells{n}.categories{c}.dFoT);
            thissigma = data.cells{n}.categories{c}.baselineSD(tr);
            thistrace = data.cells{n}.categories{c}.dFoT{tr};
            thistrace = remove_fluctuations(thistrace,75,8);
            thistrace = thistrace-median(thistrace);
            % Generate mask by transientmask()
            thismask=transientmask(thistrace,mindur,nsigma,0.5,thissigma); 
            data.cells{n}.categories{c}.transientmask{tr}=...
                reshape(thismask,[1,length(thismask)]); %170525 for consistency
        end
    end
    fprintf('Event detection for cell %1$i completed\n',n);
end
% Recalculate dFoY for significant events only.
data=recalculate_dFoY(data,true,true,0.1:0.05:2.1);
end

function mask = create_transient_mask(trace, framerate, nsigma,...
    mindur, maxdur, cutoff, GaussFiltFreq, driftdur)
% Take fluorescent trace and extract the > N sigma events at a duration that
% fullfills the requirement of less than 5% negative going events.
% Procedere: Use cat. 

% Gaussian filtering if required
if GaussFiltFreq
    gaussigma = round(framerate/GaussFiltFreq); % Create gauss filter for pre-processing
    gaussfilt = fspecial('gaussian', [10, 1], gaussigma);
else
    gaussfilt = 1;
end

trace = conv(trace, gaussfilt, 'same');

% Use Dombeck method for subtracting slow drift components
trace = subtract_drift(trace, round(driftdur*framerate), 8);

% Increase transient duration to obtain < 5% false positives
minsamples = round(mindur*framerate);
maxsamples = round(maxdur*framerate);

SD = findsigma2(trace,3); % Baseline SD excluding 3sigma events
return % DEBUG!

for len = minsamples:maxsamples
end 


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
figure; plot(bgtrace);
trace = trace-bgtrace';
figure; plot(trace-bgtrace');
end