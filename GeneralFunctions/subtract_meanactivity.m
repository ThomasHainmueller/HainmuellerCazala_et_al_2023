function data = subtract_meanactivity(data,categories,cells)
% Subtract the mean activity of all cells from each individual one. Usefull
% when there is a dominant source of signal in all data (either artefact or
% e.g. running-induced increases in fluorescence) that you want to cancel
% out. You can specify a subset of cells from which the average is to be
% made (e.g. when one or more 'background' ROIs were created).

if nargin < 3
    cells = 1:length(data.cells);
end

if nargin < 2
    categories = 1:length(data.metadata.categories);
end

% Get the data from all cells and average over them
for c = categories
    for r = 1:length(data.metadata.categories{c}.ft)
        % Get activity from all cells for this run
        activity = [];
        for n = cells
            activity(:,n) = data.cells{n}.categories{c}.dFoT{r};
        end

        % Subtract mean over all cells from each cell
        for n = cells
            data.cells{n}.categories{c}.dFoT{r}=...
                activity(:,n) - nanmean(activity,2);
        end
        %figure; plot(mean(activity,2));
        %figure; imagesc(activity');
    end
end
end