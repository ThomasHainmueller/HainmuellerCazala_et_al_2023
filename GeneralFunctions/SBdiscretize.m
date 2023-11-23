function discrete = SBdiscretize(data, coordinates, bins, mode)
% Function to discretize a data vector according to the coordinates where
% each datapoint was acquired. Takes two equal length double vectors as inputs.
% The third vector 'bins' must be a continous monotonic growing vector of type
%[ minvalue : stepsize : maxvalue ]. 'discrete' is a vector of lenght bins
%- 1 that holds the mean of the values in data that were acquired in each
%bin.
% 
% EDIT 6/10/2023: Renamed this function from 'discretize' to resolve naming
% conflict with the built in matlab function of the same name!
if nargin < 4
    mode = 'mean'; % Can be 'mean', 'median' or 'distribution', in the latter case, 
                   % the distribution of 'data' values for each bin is 
                   % returned as a cell array. Feature added 230815
end

for n = length(bins)-1:-1:1
    if strcmp(mode,'mean')
        discrete(n,1) = nanmean(data((coordinates>bins(n))&(coordinates<=bins(n+1))));
    elseif strcmp(mode,'median')
        discrete(n,1) = nanmedian(data((coordinates>bins(n))&(coordinates<=bins(n+1))));
    elseif strcmp(mode,'distribution');
        discrete{n,1} = data((coordinates>bins(n))&(coordinates<=bins(n+1)));
    else
        warning('Mode must be ''mean'' or ''distribution'', returning empty');
        discrete = [];
        return
    end
end

end