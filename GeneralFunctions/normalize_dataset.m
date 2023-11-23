function Dataset = normalize_dataset(Dataset, Fbaseline)
% Normalize all the traces in a dataset of lineartracks.

if nargin < 2
    Fbaseline = 0.2;
end

for r = length(Dataset):-1:1
    thisdata = Dataset{r}.data;  
    % Exclude fields 1:3, i.e. frametimes, x and y coords.
    for c = size(thisdata,2):-1:4
        thisdata(:,c) = normtobaseline(thisdata(:,c),Fbaseline);
    end    
    Dataset{r}.data=thisdata;
    clear thisdata
end

end