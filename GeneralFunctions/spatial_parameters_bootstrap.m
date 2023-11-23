function data = spatial_parameters_bootstrap(data,nshuffle,bins)
% Updated version of the bootstrapping. Save the random maps with the
% dataset and compute p-value for spatial info, spatial coherence, vector
% tuning.
%
% Thomas Hainmueller, 6/2023 
if nargin<3
    bins = 0.1:0.05:2.1;
end

if nargin<2
    nshuffle=1000;
end

rng('default')

for c = 1:length(data.cells{1}.categories)
    ytrace = deal(data.metadata.categories{c}.y);
    ytrace = cat(2,ytrace{:});
    moving = deal(data.metadata.categories{c}.moving);
    moving = logical(cat(2,moving{:}));
        %TODO: Fix framerate issue, this should allready be in FPS format!
        %framerate = deal(data.metadata.categories{c}.acquisition_rate);
        %framerate = 1000/mean(cat(2,framerate{:})); % frames/second
    framerate = mean(data.metadata.categories{c}.acquisition_rate);

    
    goodidcs = moving & ytrace >= bins(1) & ytrace <= bins(end);
    ytrace = ytrace(goodidcs);
    
    % SIC! Using built-in discretize function here; Gets around the issue
    % with 'edges' in the 'hist' function, see equivalent in spatial_info.m
    yBinIDs = discretize(ytrace,bins);
%     yBinBool = repmat(yBinIDs,length(bins)-1,1)';
%     yBinBool = yBinBool == 1:length(bins)-1;
    yhist = hist(yBinIDs,1:length(bins)-1)./length(yBinIDs);
    
    shifts = randsample(length(ytrace),nshuffle,true);
    
    for n = 1:length(data.cells)
        transients=deal(data.cells{n}.categories{c}.transientmask);
        transients=cat(2,transients{:});
        
        signals=deal(data.cells{n}.categories{c}.dFoT);
        signals=cat(2,signals{:}).*transients;
        
        signals = signals(goodidcs);
        dFmean = mean(signals);
        
        % 'Real' spatial information value
        %dFoY = makeMap(signals,yBinBool);
        dFoY = SBdiscretize(signals,ytrace,bins);
        Ispatial = si_on_map(dFoY,yhist,dFmean,framerate);
        
        Iboot = NaN(nshuffle,1);
        parfor it = 1:nshuffle
            %dFoY = makeMap(circshift(signals,shifts(it)),yBinBool);
            dFoY = SBdiscretize(circshift(signals,shifts(it)),ytrace,bins);
            Iboot(it) = si_on_map(dFoY,yhist,dFmean,framerate);
        end
        
        % P-value of real spatial information
        p = length(Iboot(Iboot>Ispatial))/length(Iboot);
        
        % Write to dataset:
        data.cells{n}.categories{c}.SpatialInfoBoot=Iboot;
        data.cells{n}.spatialinfo(c)=Ispatial;
        data.cells{n}.spatial_P(c)=p;
    end
end
end

% function dFoY = makeMap(signal,yBinBool)
% % Bin data with logical indexing.
% for b = 1:size(yBinBool,2)
%     dFoY(b) = nanmean(signal(yBinBool(:,b)));
% end
% end

function si = si_on_map(dFoY,yhist,dFmean,framerate)
% Simplified computation of spatial info on pre-computed map and
% parameters. See 'spatial_info.m' for full from-scratch version.

for x = length(dFoY):-1:1
    % Spatial information, Skaggs et al., 1993
    Ispatial(x) = dFoY(x)*log2(dFoY(x)/dFmean)*yhist(x);
end

si = nansum(Ispatial)*framerate;

end