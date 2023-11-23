function data = split_dataset(data,varargin)
% Extract a sub-dataset from an existing one by selecting a subset of cells
% or categories from the old file.

args=struct('categories',1:length(data.metadata.categories),...
    'cells',1:length(data.cells));
    % will appear in the order specified in input argument.

% Overwrite default parameters if required
for pair = reshape(varargin,2,[])
    args.(pair{1}) = pair{2};
end

newdata=data;

% Copy selected cells to new dataset.
if any(args.cells~=1:length(data.cells))
    cells=args.cells;
    newdata=rmfield(newdata,'cells');
    for n=1:length(cells)
        newdata.cells{n}=data.cells{cells(n)};
    end
end

% Copy category related informations to new dataset.
if any(args.categories~=1:length(data.metadata.categories))
    cats=args.categories;
    newdata.metadata=rmfield(newdata.metadata,'categories');
    for c=1:length(cats)
        % TODO regrouping of category related information in the 
        newdata.metadata.catgeories{c}=data.metadata.categories{cats(c)};

end