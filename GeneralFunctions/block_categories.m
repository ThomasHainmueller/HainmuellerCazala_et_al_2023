function newdata = block_categories(data,categories, blocksize)
% create a new datafile in which the specified categories are split into
% blocks of 'blocksize' consecutive trials and arranged in alternating
% order. This works only for two selected categories at a time!

if nargin<3
    blocksize = 5; % Number of runs going into each block
end

% Determine the (maximum) number of blocks in each category.
nruns = [];
for c = categories
    nruns(end+1) = length(data.metadata.categories{c}.ft);
end
nblocks = ceil(max(nruns)/blocksize);

% Match category numbers between old and new ds
oldcats(1:2:2*nblocks) = categories(1);
oldcats(2:2:2*nblocks) = categories(2);
newcats = 1:2*nblocks;
%newdata = regroup_categories(data, 1:2*nblocks, oldcats);

% Create new dataset, transferring all selected ('hard coded') fields
try
    newdata.label = data.label;
end

try
    newdata.animal = data.animal;
end

% 1) Metadata, general part
newdata.metadata.Placefield_P(:,1:2*nblocks)=...
    data.metadata.Placefield_P(:,oldcats);

% 2) Cells, general part
for n = 1:length(data.cells)
    try
        newdata.cells{n}.roi = data.cells{n}.roi;
    end
    newdata.cells{n}.transient_sigma = data.cells{n}.transient_sigma;
    newdata.cells{n}.transient_mindur = data.cells{n}.transient_mindur;
    newdata.cells{n}.transientrate(newcats) = data.cells{n}.transientrate(oldcats);
    newdata.cells{n}.AUCrate(newcats) = data.cells{n}.AUCrate(oldcats);
    newdata.cells{n}.spatialinfo(newcats) = data.cells{n}.spatialinfo(oldcats);
    newdata.cells{n}.spatial_P(newcats) = data.cells{n}.spatial_P(oldcats);
    newdata.cells{n}.Placefield_P(newcats) = data.cells{n}.Placefield_P(oldcats);
end
    
for c = newcats
    % Get indices for the trial block to be deposited in this category.
    block = ceil(c/2); % For two categories
    if length(data.metadata.categories{oldcats(c)}.ft) < (block-1)*blocksize+1
        break
    elseif length(data.metadata.categories{oldcats(c)}.ft) < block*blocksize
        theseruns = (block-1)*blocksize+1 : length(data.metadata.categories{oldcats(c)}.ft);
    else
        theseruns = (block-1)*blocksize+1 : block*blocksize;
    end
    
    % 3) Metadata, categorial part
    try % Optional part
        newdata.metadata.categegories{c}.rewardszone = ...
            data.metadata.categories{oldcats(c)}.rewardzone;
        newdata.metadata.categegories{c}.rewardpoints = ...
            data.metadata.categories{oldcats(c)}.rewardpoints;
    end
    newdata.metadata.categories{c}.acquisition_rate = ...
        data.metadata.categories{oldcats(c)}.acquisition_rate(theseruns);
    newdata.metadata.categories{c}.ft = ...
        data.metadata.categories{oldcats(c)}.ft(theseruns);    
    newdata.metadata.categories{c}.x = ...
        data.metadata.categories{oldcats(c)}.x(theseruns);    
    newdata.metadata.categories{c}.y = ...
        data.metadata.categories{oldcats(c)}.y(theseruns);
    try
        newdata.metadata.categories{c}.licks = ...
            data.metadata.categories{oldcats(c)}.licks(theseruns);
        newdata.metadata.categories{c}.pupil_area = ...
            data.metadata.categories{oldcats(c)}.pupil_area(theseruns);
    end
    try
        newdata.metadata.categories{c}.filename = ...
            data.metadata.categories{oldcats(c)}.filename(theseruns);
    end
    newdata.metadata.categories{c}.moving = ...
        data.metadata.categories{oldcats(c)}.moving(theseruns);    
    try
        newdata.metadata.categories{c}.tforms = ...
            data.metadata.categories{oldcats(c)}.tforms(theseruns);
    end
    
    % Cell data, categorial part:
    for n = 1:length(data.cells)
        newdata.cells{n}.categories{c}.dFoT = ...
            data.cells{n}.categories{oldcats(c)}.dFoT(theseruns);
        try
            newdata.cells{n}.categories{c}.F0 = ...
                data.cells{n}.categories{oldcats(c)}.F0(theseruns);
        end
        newdata.cells{n}.categories{c}.baselineSD = ...
            data.cells{n}.categories{oldcats(c)}.baselineSD(theseruns);
        newdata.cells{n}.categories{c}.dFoY = ...
            data.cells{n}.categories{oldcats(c)}.dFoY(theseruns);
        newdata.cells{n}.categories{c}.transientmask = ...
            data.cells{n}.categories{oldcats(c)}.transientmask(theseruns);
        newdata.cells{n}.categories{c}.SpatialInfoBoot = ...
            data.cells{n}.categories{oldcats(c)}.SpatialInfoBoot;
        newdata.cells{n}.categories{c}.placefield = ...
            data.cells{n}.categories{oldcats(c)}.placefield;
    end
end
end