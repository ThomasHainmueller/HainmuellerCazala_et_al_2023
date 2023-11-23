function data = PVcorrelation(data,cells)
% Append a matrix with the population vectors (for each category) and the
% population vector correlations to the dataset (metadata section). Use the
% dFoY traces stored with the file allready.
if nargin<2
    cells = 1:length(data.cells); % Default: Calculate corr for all cells
end

for c = length(data.cells{1}.categories):-1:1
    for n = length(data.cells):-1:1;
        for tr = 1:length(data.cells{n}.categories{c}.dFoY)
            thisdFoY(tr,:) = data.cells{n}.categories{c}.dFoY{tr};
        end
        PV(n,c,:) = nanmean(thisdFoY,1);
    end
    data.metadata.categories{c}.PV = squeeze(PV(:,c,:)); % Store individulay w. category
    %figure; imagesc(squeeze(PV(:,c,:)));
end

for bin = size(PV,3):-1:1
    %[r(:,:,n),p]=corrcoef(dFoY,'rows','complete');
    r(:,:,bin) = corrcoef(PV(cells,:,bin),'rows','complete');
end

data.metadata.PVcorr = nanmean(r,3); % cxc matrix with PV correlations
%figure; imagesc(data.metadata.PVcorr); colormap('jet'); caxis([0 1]);
end
            
