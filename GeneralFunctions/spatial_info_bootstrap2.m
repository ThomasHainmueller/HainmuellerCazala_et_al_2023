function data = spatial_info_bootstrap2(data,nshuffle,ywindow,bins)
% Updated version with faster core algorithm and multi-threading 180320
%
% Take a dataset and calculate for each category and each cell the spatial
% information represented as well as the distribution of spatial
% informations with shuffled y traces. Return spatial info (in bits/s)
% and p-value for the bootstraped distribution.
% By default only for moving periods, can be easily adapted (see
% transientrates() and so on).

if nargin<4
    bins = 0.1:0.05:2.1;
end
if nargin<3
    ywindow=50; % 50 Frames / approx. 10 s
end
if nargin<2
    nshuffle=1000;
end

for c = 1:length(data.cells{1}.categories)
    % Rationale: Spatial info derives from mean firing rate per bin and
    % probability of occupancy (p). p is the same for all cells and for all
    % permutations of the spatial trajectory. The permutations are the same
    % for all cells. To expedite, the permuted y-traces are replaced with
    % lookup-tables (logical of dim(nsamples x nbins)).
    
    ytrace = deal(data.metadata.categories{c}.y);
    ytrace = cat(2,ytrace{:});
    moving = deal(data.metadata.categories{c}.moving);
    moving = logical(cat(2,moving{:}));
    framerate=mean(data.metadata.categories{c}.acquisition_rate);
    
    % bootstrapped y traces (moving periods only!)
    yrand = shuffle(ytrace(moving),nshuffle,ywindow);
    
    for b = length(bins)-1:-1:1
        % Calculate probability of mouse occupying each bin
        occupancy(b) = length(find(...
            (ytrace>bins(b)) & (ytrace<=bins(b+1)))) / length(ytrace);
        for rd = nshuffle:-1:1
            % Order 1=nshuffles, 2=bins, 3=frames
            ylookup(rd, b, :) = (yrand(:,rd)>bins(b)) & (yrand(:,rd)<=bins(b+1));
        end
    end
    
    for n = 1:length(data.cells)
        transients=deal(data.cells{n}.categories{c}.transientmask);
        transients=cat(2,transients{:});
        signals=deal(data.cells{n}.categories{c}.dFoT);
        signals=cat(2,signals{:}).*transients;
        signals = signals(moving);
        meanrate = nanmean(signals);
        
        % 'Real' spatial information value
        Ispatial = spatial_info(signals, ytrace(moving), bins, framerate);
        
        % spatial info for bootstraped ytraces
         %for ytr = nshuffle:-1:1 
             %Iboot(ytr)=spatial_info(signals(moving),yrand(:,ytr),bins,framerate);
         %end
        for rd = nshuffle:-1:1
            for b = length(bins)-1:-1:1
            	thisrate = nanmean(signals(ylookup(rd, b, :)));
                
                % Spatial information, Skaggs et al., 1993
                thisIspatial(b) = thisrate*log2(thisrate/meanrate)*occupancy(b);
            end
            Iboot(rd) = nansum(thisIspatial)/length(signals)*framerate;
        end
        
        % P-value of real spatial information
        p = length(Iboot(Iboot>Ispatial))/length(Iboot);
        
        % Write to dataset:
        data.cells{n}.categories{c}.SpatialInfoBoot=Iboot;
        data.cells{n}.spatialinfo(c)=Ispatial;
        data.cells{n}.spatial_P(c)=p;
    end
    clear yrand ylookup
end
end