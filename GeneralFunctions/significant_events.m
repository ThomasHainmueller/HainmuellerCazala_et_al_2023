function data = significant_events(data, sigmarange, timerange, cutoff)
% Takes a dataset created by get_lineartracks -> convert_data, defines the
% minimum elevation (in sigma of baseline) and duration (in frames) that
% defines a significant transient (Ratio between positive and negative Transients
% <0.05; see Dombeck et al., 2007). Creates masks for all traces and stores
% them with the dataset individually for all cells in the dataset.
if nargin < 4
    cutoff = 0.05;
end
if nargin < 3
    timerange = 1:5;
end
if nargin < 2
    sigmarange = 2:0.5:5;
end

% Loop over cells in dataset
for n = 1:length(data.cells)
    % Collect and collate all calcium transients from this cell
    thisdFoT = [];
    for c = 1:length(data.cells{n}.categories)
        traces = deal(data.cells{n}.categories{c}.dFoT);
        traces = cat(2,traces{:});
        % Subtract the 8% value for each 75 frames (usually ~15s)
        traces = remove_fluctuations(traces,75,8);
        % Subtract median to assure that baseline is centered to zero
        traces = traces-median(traces);
        thisdFoT = cat(2,thisdFoT,traces);
        % Note the standard deviation of baseline for each category
        sigma(c)=findsigma(traces,100,10);
        clear traces
    end
    
    % Find the values for nsigma and dT that yield the highest detection
    % rate within the cutoff tolerance of error rate.
    [nsigma, dT, pval]=significant_event_par(thisdFoT,sigmarange,timerange,cutoff);
    data.cells{n}.transient_sigma = nsigma;
    data.cells{n}.transient_dT = dT;
    data.cells{n}.transient_pval = pval;
    
    % Create masks for all traces in all categories and store them with the
    % dataset.
    for c = 1:length(data.cells{n}.categories)
        data.cells{n}.categories{c}.baselinesigma = sigma(c);
        for tr = 1:length(data.cells{n}.categories{c}.dFoT);
            thistrace = data.cells{n}.categories{c}.dFoT{tr};
            thistrace = remove_fluctuations(thistrace,75,8);
            thistrace = thistrace-median(thistrace);
            % Generate mask by transientmask()
            thismask=transientmask(thistrace,dT,nsigma,0.5,sigma(c)); 
            data.cells{n}.categories{c}.transientmask{tr}=...
                reshape(thismask,[1,length(thismask)]); %170525 dimension consitency
        end
    end
    fprintf('Event detection for cell %1$i completed\n',n);

end
% Recalculate dFoY for significant events only.
%data=recalculate_dFoY(data,true,true,0.1:0.025:2.1);
end