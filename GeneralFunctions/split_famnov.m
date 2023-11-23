function newdata = split_famnov(data,varargin)
% utility function for interneuron datasets where familiar and novel
% contexts were acquired in the same recording files (rapid switching).
% Splits traces in the respective recordings into the familiar and novel
% segment as indicated by the x-coordinate and stores them into separate
% categories of the dataset. All other categories are incremented by one.

args=struct('category',1,'splitlevel',0); 
% splitlevel: by default (0), the half-maximum x-value of the first dataset 
% is used as category divider. May be set manually if required.

% Overwrite default parameters if required
for pair = reshape(varargin,2,[])
    %args.(pair{1}) = pair{2};
    if isfield(args,pair{1})
        args.(pair{1})=pair{2};
    else
        error('Input argument does not exist\n')
    end
end

% Create new dataset with room for the additional split categories.
% Special case: If there is just 1 category, then duplicate it.
if ~(length(data.metadata.categories)==1)
    newdata = regroup_categories(data, 1:length(data.metadata.categories)+1,...
        [1:args.category,0,args.category+1:length(data.metadata.categories)]);
else
    newdata = regroup_categories(data,[1,2],[1,1]);
end

% Infer level (x-coordinate) at which to split the familiar and novel part.
if ~args.splitlevel
    args.splitlevel = max(data.metadata.categories{args.category}.x{1})/2;
end
    
% Iterate over traces and split each one in familiar and novel segment.
for tr = 1:length(data.metadata.categories{args.category}.x)
    
    % Split all metadata traces into the familiar and novel segment. 
    % Default filtering at 50 frames to erase artificial fluctuations
    thisx = medfilt1(single(data.metadata.categories{args.category}.x{tr}), 100); 
    divider = thisx>=args.splitlevel;
    
    newdata.metadata.categories{args.category}.ft{tr}=...
        data.metadata.categories{args.category}.ft{tr}(~divider);
    newdata.metadata.categories{args.category+1}.ft{tr}=...
        data.metadata.categories{args.category}.ft{tr}(divider);  
    newdata.metadata.categories{args.category}.x{tr}=...
        data.metadata.categories{args.category}.x{tr}(~divider);
    newdata.metadata.categories{args.category+1}.x{tr}=...
        data.metadata.categories{args.category}.x{tr}(divider);
    newdata.metadata.categories{args.category}.y{tr}=...
        data.metadata.categories{args.category}.y{tr}(~divider);
    newdata.metadata.categories{args.category+1}.y{tr}=...
        data.metadata.categories{args.category}.y{tr}(divider);
    newdata.metadata.categories{args.category}.moving{tr}=...
        data.metadata.categories{args.category}.moving{tr}(~divider);
    newdata.metadata.categories{args.category+1}.moving{tr}=...
        data.metadata.categories{args.category}.moving{tr}(divider);
    
    % Split cellular traces into the familiar and novel segment
    for n = 1:length(data.cells)
        newdata.cells{n}.categories{args.category}.dFoT{tr}=...
            data.cells{n}.categories{args.category}.dFoT{tr}(~divider);
        newdata.cells{n}.categories{args.category+1}.dFoT{tr}=...
            data.cells{n}.categories{args.category}.dFoT{tr}(divider);
    end
end

% Propagate other parameters to the novel category
newdata.metadata.categories{args.category+1}.acquisition_rate=...
    data.metadata.categories{args.category}.acquisition_rate;
newdata.metadata.categories{args.category+1}.filename=...
    data.metadata.categories{args.category}.filename;

for n = 1:length(data.cells)
    newdata.cells{n}.categories{args.category+1}.F0=...
        data.cells{n}.categories{args.category}.F0;
    newdata.cells{n}.categories{args.category+1}.baselineSD=...
        data.cells{n}.categories{args.category}.baselineSD;
end

newdata = recalculate_dFoY(newdata,0,0);

end

