function data = subtract_background(data,bg)
% Subtract background from interneuron datasets

allfn = {};
%fnIDCs = [];
for c = 1:length(data.metadata.categories)
    allfn = [allfn, data.metadata.categories{c}.filename];
%     theseIDCs = [ones(1,length(data.metadata.categories{c}.filename))*c;...
%         1:length(data.metadata.categories{c}.filename)];
%     fnIDCs = cat(2,fnIDCs,theseIDCs);
end

if length(allfn) ~= length(bg.filename)
    warning('Background not available for all files in session, aborting...');
    return
end

for c = 1:length(data.metadata.categories)
    for r = 1:length(data.metadata.categories{c}.filename)
        % Find correct background
        bgID = find(strcmp([data.metadata.categories{c}.filename{r} '.sbx'],bg.filename));
    
        % Subsample background to match data
        xq = linspace(1,length(bg.trace{bgID}),...
            length(data.metadata.categories{c}.ft{r}));
        thisbg = interp1(1:length(bg.trace{bgID}),bg.trace{bgID},xq);
        
        for n = 1:length(data.cells)
            nbg = thisbg / data.cells{n}.categories{c}.F0{r};
            data.cells{n}.categories{c}.dFoT{r} = data.cells{n}.categories{c}.dFoT{r} -...
                nbg + mean(nbg);
        end
    end
end
          
end