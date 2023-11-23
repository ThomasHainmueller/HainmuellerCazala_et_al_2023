function plot_placemaps(data,indices,sortcategory,dFoT,transients,active,SI,PF,zscored)
% Plots lineartrack acquired activity grouped to category. As a summary for
% all cells in the dataset. dFoT: Toggle all traces over time (1/0)
if nargin < 9
    zscored = false;
end
if nargin < 8
    % Place field significance from bootstrapping
    PF = 0; % 0.05 = significant Placefields only
end
if nargin < 7
    % Cells whos spatial information exceeds the given significance.
    SI = 0; % 0=Criterion does not apply, 0.05 is standard.
end
if nargin < 6
    % Activity level. Default = 1/60 transients (i.e. 1/min).
    active = 0; % 0=Criterion not applied
end
if nargin < 5
    % Use only significant transients
    transients = false;
end
if nargin < 4
    % Plot all traces in the bottom plot over time
    dFoT = true;
end

if nargin < 3
    sortcategory = 1;
end

if nargin < 2
    % Default: plot all cells in dataset.
    indices = 1:length(data.cells);
end

% Convert indices to logical array
Lindices=false(1,length(data.cells));
Lindices(indices)=true;

% Find active and spatially modulated cells
if active
    activecells=false(1,length(data.cells));
    for n =1:length(data.cells)
        if strcmp(class(data.cells{n}.transientrate),'double') % Downward compatibility
            if data.cells{n}.transientrate(sortcategory)>=active
                activecells(n)=true;
            end
        else
            if data.cells{n}.transientrate{sortcategory}>=active
                activecells(n)=true;
            end
        end
    end
    Lindices = Lindices&activecells;
end
if SI
    SIcells=false(1,length(data.cells));
    for n=1:length(data.cells)
        if SI == 1
            SIcells(n)=true; % Exclude criterion
        elseif strcmp(class(data.cells{n}.transientrate),'double') % Downward compatibility
            if data.cells{n}.spatial_P(sortcategory)<=SI
                SIcells(n)=true;
            end
        else
            if data.cells{n}.spatial_P{sortcategory}<=SI
                SIcells(n)=true;
            end
        end
    end
    Lindices = Lindices&SIcells;
end
if PF
    PFcells=false(1,length(data.cells));
    for n=1:length(data.cells)
        if data.cells{n}.Placefield_P(sortcategory)<=PF
            PFcells(n)=true;
        end
    end
    Lindices = Lindices&PFcells;
end

filter = fspecial('gaussian',[9,1],1.0);

for c = length(data.metadata.categories):-1:1
    cats{c}.alltracks = [];
    newcellno = 1;
    % Loop over cells and collect activities over time and Y
    for cell = find(Lindices)
        thisdFoT = [];
        thisdFoY = [];
        %thismask = [];
        for run = 1:length(data.cells{cell}.categories{c}.dFoY) 
            if transients
                thismask = reshape(data.cells{cell}.categories{c}.transientmask{run},...
                    size(data.cells{cell}.categories{c}.dFoT{run}'));
                %thismask = data.cells{cell}.categories{c}.transientmask{run}; %Removed matrix inversion 170519
                if zscored
                    thisdFoT = cat(1,thisdFoT,data.cells{cell}.categories{c}.zscored{run}'.*thismask);
                else
                    thisdFoT = cat(1,thisdFoT,data.cells{cell}.categories{c}.dFoT{run}'.*thismask);
                end
                rawdFoY = data.cells{cell}.categories{c}.dFoY{run};
            else
                if zscored
                    thisdFoT = cat(1,thisdFoT,data.cells{cell}.categories{c}.zscored{run}');
                    rawdFoY = SBdiscretize(data.cells{cell}.categories{c}.zscored{run},...
                        data.metadata.categories{c}.y{run},0.1:0.025:2.1);
                else
                    thisdFoT = cat(1,thisdFoT,data.cells{cell}.categories{c}.dFoT{run}');
                    rawdFoY = SBdiscretize(data.cells{cell}.categories{c}.dFoT{run},...
                        data.metadata.categories{c}.y{run},0.1:0.025:2.1);
                end
            end
            % Use Gaussian filter to smooth the space plot.
            nanvals = isnan(rawdFoY);
            rawdFoY(nanvals)=0;
            rawdFoY = conv(rawdFoY,filter,'same');
            rawdFoY(nanvals)=NaN;
            thisdFoY = cat(2,thisdFoY,rawdFoY);
        end
        % Append dFoT activity from this cell for all trials.
        cats{c}.alltracks = cat(2,cats{c}.alltracks,thisdFoT);
        
        % Average dFoY activity over all trials
        cats{c}.spaceplot(:,newcellno) = nanmean(thisdFoY,2);
        newcellno=newcellno+1;
        clear thisdFoT thisdFoY
    end
end

% Find order of placefileds in selected category
[~, sortvector] = max(cats{sortcategory}.spaceplot);
[~, sortindices] = sort(sortvector);


for c = 1:length(data.metadata.categories)
    if dFoT
    FigHandle = figure('Position',[(c-1)*350,50,500,900]);
    subplot(2,1,1);
        imagesc(cats{c}.spaceplot(:,sortindices)');
        title(['Category ' num2str(c)])
        colormap('jet');
        if zscored
            caxis([0.0 2.5]);
        else
            caxis([0.0 0.3]);
        end
    subplot(2,1,2);
        imagesc(cats{c}.alltracks(:,sortindices)');
        colormap('jet');
        if zscored
            caxis([0.0 8.0]);
        else
            caxis([0.0 2.0]);
        end
    else
    FigHandle = figure('Position',[(c-1)*450,50,500,900]);
    imagesc(cats{c}.spaceplot(:,sortindices)');
    title(['Category ' num2str(c)])
    caxis([0.0 0.3]);
    colormap('jet');
    end
end
end
