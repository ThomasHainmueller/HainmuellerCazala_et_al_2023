function [cats sortindices]= plot_allcells(dataset,indices, sortcategory,dFoT)
% Plots lineartrack acquired activity grouped to category. As a summary for
% all cells in the dataset. dFoT: Toggle all traces over time (1/0)

if nargin < 4
    dFoT = true;
end

if nargin < 3
    sortcategory = dataset{1}.category;
end

if nargin < 2
    % Default: plot all cells in dataset.
    indices = [4:length(dataset{1}.data(:,1))];
end

bins = [0.05:0.025:2.1];
numruns = length(dataset);
numcells = length(indices);
%numcategories = 12;

%res = double.empty(length(bins)-1,numcells,0,numruns);

% Assign the indices of all runs to their respective group.
for run = 1:numruns
    thiscat = dataset{run}.category;
    try
        cats{thiscat}.runs = [cats{thiscat}.runs, run];
    catch
        cats{thiscat}.runs = [run];
    end
end

% remove emtpty categories.
cats = cats(~cellfun('isempty',cats));
numcategories = length(cats);

for c = numcategories:-1:1
    cats{c}.alltracks = [];
    cats{c}.spaceplot = double.empty(length(bins)-1,numcells,0);
    
    for r = cats{c}.runs        
        % Make contionuous plot of traces over time
        cats{c}.alltracks = cat(2,cats{c}.alltracks,dataset{r}.data);
        thisy=dataset{r}.data(3,:);
        Nrun=size(cats{c}.spaceplot,3)+1;
        
        % make dF/distance on track plots for all cells
        for cellno = 1:length(indices)
            oldID = indices(cellno);
            % CAVE: dF/F baseline set to 0 here!
            % EDITED 15/12/05 
            thistrace=dataset{r}.data(oldID,:);
            %thistrace = normtobaseline(dataset{r}.data(oldID,:),0.1);
            cats{c}.spaceplot(:,cellno,Nrun) = SBdiscretize(thistrace,thisy,bins);
            
            clear thistrace
        end
    end
    % Average placeactivity over all trials
    cats{c}.spaceplot(cats{c}.spaceplot==0) = NaN;
    cats{c}.spaceplot = nanmean(cats{c}.spaceplot,3);
end

% Find order of placefields in selected category
[~, sortvector] = max(cats{sortcategory}.spaceplot);
[~, sortindices] = sort(sortvector);

for c = 1:numcategories
    if dFoT
    FigHandle = figure('Position',[(c-1)*450,50,500,900]);
    subplot(2,1,1);
        imagesc(cats{c}.spaceplot(:,sortindices)');
        title(['Category ' num2str(c)])
        caxis([1.0 1.5]);
    subplot(2,1,2);
        imagesc(cats{c}.alltracks(sortindices+3,:));
        caxis([1.0 3.0]);
    else
    FigHandle = figure('Position',[(c-1)*450,50,500,900]);
    imagesc(cats{c}.spaceplot(:,sortindices)');
    title(['Category ' num2str(c)])
    caxis([1.0 1.5]);
    end
end
end
