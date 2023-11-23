function dataset = convert_data(originaldata, bins, minspeed, PVcells)
% Convert a dataset created by get_lineartracks to a cell by cell
% representation of calcium signals with the following structure:
% data.cells{1:n}.categories{1:m}.runs{1:o}.dFoT,dFoY
% data.metadata.categories{1:m}.runs{1:o}.x,y,ft,movingmask
%
% Parameters:
% bins = bin size for dFoY discretization
% minspeed = 4m/5Hz: 5E-3; 4m/15Hz: 1.67E-3; 4m/30Hz: 8.33E-4; 
%       minimum (average) speed for excluding non moving periods
% PVcells = special conversion for tdT/GCaMP double labeled cells

if nargin < 4
    PVcells = false;
end
if nargin < 3
    minspeed = 1.67E-3; % Approx. 5cm/sec at V2.1=4m and 5Hz framerate
end
if nargin < 2
    % bin size = 5cm for a 4m virtual track.
    bins = 0.05:0.025:2.1;
end

ncells = size(originaldata{1}.data,1)-3;

% Distribute traces from the original dataset upon the new structure
for run = 1:length(originaldata)
    thiscat = originaldata{run}.category;
    
    % find out how many runs have been stored in this category.
    try
        catrun = length(dataset.metadata.categories{thiscat}.ft)+1;
    catch
        catrun = 1;
    end
    
    dataset.metadata.categories{thiscat}.ft{catrun} = originaldata{run}.data(1,:);
    dataset.metadata.categories{thiscat}.x{catrun} = originaldata{run}.data(2,:);
    dataset.metadata.categories{thiscat}.y{catrun} = originaldata{run}.data(3,:);
    dataset.metadata.categories{thiscat}.filename{catrun}=...
        originaldata{run}.filename;
    dataset.metadata.categories{thiscat}.moving{catrun}=...
        movingmask(originaldata{run}.data(3,:),minspeed);
    dataset.metadata.categories{thiscat}.acquisition_rate(catrun)=...
        1000/mean(diff(originaldata{run}.data(1,:))); % In Hz
    
    % Append background images and ROIs if existant
    if isfield(originaldata{run},'greenavg')
        dataset.metadata.categories{thiscat}.greenavg{catrun}=...
            originaldata{run}.greenavg;
    end
    if isfield(originaldata{run},'redavg')
        dataset.metadata.categories{thiscat}.redavg{catrun}=...
            originaldata{run}.redavg;
    end
    if isfield(originaldata{run},'roimask')
        dataset.metadata.categories{thiscat}.roimask{catrun}=...
            originaldata{run}.roimask;
    end 
    
    for cell = ncells:-1:1
        % Calculate baseline and bl sigma of the trace without transients
        if PVcells
            thisdFoT=originaldata{run}.data(cell+3,:)-1.0;
            thisGraw=originaldata{run}.greenraw(cell+3,:);
            thisRraw=originaldata{run}.redraw(cell+3,:);
            dataset.cells{cell}.categories{thiscat}.greenraw{catrun}=thisGraw;
            dataset.cells{cell}.categories{thiscat}.redraw{catrun}=thisRraw;
            dataset.cells{cell}.categories{thiscat}.meanGoR(catrun)=...
                mean(thisGraw)/mean(thisRraw);
        else  
            [thisdFoT,F0,sigma] = normtobaseline2(originaldata{run}.data(cell+3,:),2,5);
            dataset.cells{cell}.categories{thiscat}.baseline(catrun) = F0;
            dataset.cells{cell}.categories{thiscat}.baselineSD(catrun) = sigma;
        end
        thismask = dataset.metadata.categories{thiscat}.moving{catrun};
        dataset.cells{cell}.categories{thiscat}.dFoT{catrun} = thisdFoT;
        dataset.cells{cell}.categories{thiscat}.dFoY{catrun} =...
            SBdiscretize(thisdFoT.*thismask,originaldata{run}.data(3,:),bins);
        % Calculate dFodY (running speed), default 0-32cm/s range is
        % assumed; For the entire (non-trace based) Calcium signal.
        dataset.cells{cell}.categories{thiscat}.dFodY{catrun} =...
            dFoSpeed(thisdFoT,originaldata{run}.data(3,:),0:1.25E-3:4E-2);
    end
end
end