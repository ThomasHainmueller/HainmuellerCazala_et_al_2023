function data = correct_time(data)
% In some datasets, length of ft, x, y and dFoT seems to vary by one
% datapoint. Iterate through the dataset and enforce compliance by cutting
% to the length of the shortest vector.

metafields = {'ft','x','y','moving'};
cellfields = {'dFoT','transientmask'};

for c = 1:length(data.metadata.categories)
    for r = 1:length(data.metadata.categories{c}.ft)
        l = [];
        for mf = 1:length(metafields)
            l = [l, length(data.metadata.categories{c}.(metafields{mf}){r})];
        end
        for cf = 1:length(cellfields)
            l = [l, length(data.cells{1}.categories{c}.(cellfields{cf}){r})];
        end
        l = min(l);
        
        % Conversion
        for mf = 1:length(metafields)
            data.metadata.categories{c}.(metafields{mf}){r} = ...
                data.metadata.categories{c}.(metafields{mf}){r}(1:l);
        end
        for cf = 1:length(cellfields)
            for n = 1:length(data.cells)
                data.cells{n}.categories{c}.(cellfields{cf}){r} = ...
                    data.cells{n}.categories{c}.(cellfields{cf}){r}(1:l);
            end
        end
    end
end