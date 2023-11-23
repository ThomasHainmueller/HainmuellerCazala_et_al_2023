function data = crop_traces(data,category,keep)
% Crop the traces of a dataset for a given category to the intervals
% specified in 'keep'. Usefull e.g. for removing blank periods, such as the
% switch from familiar to novel in the interneuron experiments.
% 
% arguments:
% keep: 2 x nruns matrix specifying the start- and endpoint of the new
% trace

if nargin<3
% manual selection of 'keep' periods
    for r = 1:length(data.metadata.categories{category}.ft)
        figure; hold on;
        plot(data.metadata.categories{category}.x{r});
        plot(data.metadata.categories{category}.y{r});
        keep(1,r) = input('Start: ');
        keep(2,r) = input('End: ');
        close all;
    end
end

for r = 1:length(data.metadata.categories{category}.ft)
    % Crop metadata
    data.metadata.categories{category}.ft{r}=...
        data.metadata.categories{category}.ft{r}(keep(1,r):keep(2,r));
    data.metadata.categories{category}.x{r}=...
        data.metadata.categories{category}.x{r}(keep(1,r):keep(2,r));
    data.metadata.categories{category}.y{r}=...
        data.metadata.categories{category}.y{r}(keep(1,r):keep(2,r));
    data.metadata.categories{category}.moving{r}=...
        data.metadata.categories{category}.moving{r}(keep(1,r):keep(2,r));
    
    % Crop cell traces
    for n=1:length(data.cells)
        data.cells{n}.categories{category}.dFoT{r}=...
            data.cells{n}.categories{category}.dFoT{r}(keep(1,r):keep(2,r));
    end 
end
end