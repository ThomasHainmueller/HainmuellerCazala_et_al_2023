function data = random_spatial_corr(data,sigma,active,pSI,pPF)
% Calculate the pairwise correlation of place maps between each cell
% and a randomly selected one from the same dataset for statistics. 
% Uses the dFoY traces currently stored with the file. Append a
% categories^2 matrix to each cell{} with the correlations and the
% corresponding p-values. Optional Gaussian smoothing with sigma.

if nargin<5
    pPF=0.05; % changed 170929
end
if nargin<4
    %pSI=0.05; % significance threshold for spatial modulation.
    pSI = 1; %changed 170929 now selected by PF rather
end
if nargin<3
    %active=1/60; % minimum activity; default 1/min
    active = 0; % take plcs independent of activitiy
end
if nargin<2
    sigma=1; % Gaussian function with sigma = 1 frame
end

if sigma
    filter = fspecial('gaussian',[6*sigma,1],sigma);
end

% Find the correlations between categories for each cell individually.

for n=length(data.cells):-1:1
    
    % Select a random set of cells for comparison in each category
    randset = randsample(1:length(data.cells), length(data.cells{n}.categories));
    
    for c=length(data.cells{n}.categories):-1:1
        % Get mean dFoY trace from a random cell in this category
        thisdFoY=deal(data.cells{randset(c)}.categories{c}.dFoY);
        thisdFoY=cat(2,thisdFoY{:});
        thisdFoY=nanmean(thisdFoY,2);

        if sigma
            nanvals = isnan(thisdFoY);
            thisdFoY(nanvals)=0;
            thisdFoY=conv(thisdFoY,filter,'same');
            thisdFoY(nanvals)=NaN;
        end
        dFoY(:,c)=thisdFoY;
    end
    [r(:,:,n),p]=corrcoef(dFoY,'rows','complete');
    data.cells{n}.Random_Pearson_r = r(:,:,n);
    data.cells{n}.Random_Pearson_p = p;
end

% Calculate mean correlations for the entire dataset.
data.metadata.SpatialCorr_all=nanmean(r,3);
if isfield(data.metadata,'Random_SpatialCorr_sign')
    data.metadata = rmfield(data.metadata,'Random_SpatialCorr_sign'); % Remove outdated results
end
if isfield(data.metadata,'Random_SpatialCorr_n')
    data.metadata = rmfield(data.metadata,'Random_SpatialCorr_n'); % Remove outdated results
end
for c = length(data.metadata.categories):-1:1
    sign = findcells(data,c,active,pSI,pPF);
    data.metadata.Random_SpatialCorr_sign(:,c) = nanmean(r(:,c,sign),3);
    data.metadata.Random_SpatialCorr_n(c) = length(sign);
end
end