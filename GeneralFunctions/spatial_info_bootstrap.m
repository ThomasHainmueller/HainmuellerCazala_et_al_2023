function data = spatial_info_bootstrap(data,nshuffle,ywindow,bins)
% Take a dataset and calculate for each category and each cell the spatial
% information represented as well as the distribution of spatial
% informations with shuffled y traces. Return spatial info (in %dFoF AUC/s)
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

rng('default');

for c = 1:length(data.cells{1}.categories)
    ytrace = deal(data.metadata.categories{c}.y);
    ytrace = cat(2,ytrace{:});
    moving = deal(data.metadata.categories{c}.moving);
    moving = logical(cat(2,moving{:}));
        %TODO: Fix framerate issue, this should allready be in FPS format!
        %framerate = deal(data.metadata.categories{c}.acquisition_rate);
        %framerate = 1000/mean(cat(2,framerate{:})); % frames/second
    framerate=mean(data.metadata.categories{c}.acquisition_rate);
    % bootstrapped y traces.
    yrand = shuffle(ytrace(moving),nshuffle,ywindow);
   
    for n = 1:length(data.cells)
        transients=deal(data.cells{n}.categories{c}.transientmask);
        transients=cat(2,transients{:});
        signals=deal(data.cells{n}.categories{c}.dFoT);
        signals=cat(2,signals{:}).*transients;
        
        % 'Real' spatial information value
        Ispatial = spatial_info(signals(moving),ytrace(moving),bins,framerate);
        
        Iboot = NaN(nshuffle,1);
        % spatial info for bootstraped ytraces
        parfor ytr = 1:nshuffle
            Iboot(ytr)=spatial_info(signals(moving),yrand(:,ytr),bins,framerate);
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