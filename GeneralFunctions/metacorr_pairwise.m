function [res,N] = metacorr_pairwise(dataset, varargin)
% Get correlations of spatial maps for cells between categories. Include
% cells that fulfill crieteria in both ('and') or either ('or') category to
% be compared. Note, in contrast to metacorr, this function yields a
% symmetric correlation matrix. N will give the number of valid (not NaN)
% observations.

p = inputParser;
addParameter(p,'parameter','Pearson_r',@ischar); % can be 'Pearson_r' or 'Random_Pearson_r'
addParameter(p,'categories',1:length(dataset{1}.metadata.categories),@isnumeric);
addParameter(p,'active',0,@isnumeric);
addParameter(p,'SI',1,@isnumeric); % minimum p-value for spatial info
addParameter(p,'PF',1,@isnumeric); % minimum p-value for place-field
addParameter(p,'region','',@ischar);
addParameter(p,'combination','or',@ischar);
addParameter(p,'plotresults',true,@islogical);

parse(p,varargin{:})

parameter = p.Results.parameter;
categories = p.Results.categories;
active = p.Results.active;
SI = p.Results.SI;
PF = p.Results.PF;
region = p.Results.region;
combination = p.Results.combination;
plotresults = p.Results.plotresults;

% Get pairwise correlations between categories for cells meeting criteria
% in either or both:

for r = 1:length(categories)
    for c = 1:length(categories)
        thisCorr = metastatistics(dataset,'cell',parameter,'active',active,...
            'SI',SI,'PF',PF,'categories',[r c],'rescategories',[r c],'region',region,...
            'combination',combination,'plotting',false);
        if ~isempty(thisCorr)
            resC{r,c} = thisCorr(:,1,2);
            N(r,c) = length(find(~isnan(thisCorr(:,1,2))));
        else
            resC{r,c} = [];
            N(r,c) = 0;
        end
    end
end

res = catuneven(resC(:),NaN)';
res = reshape(res,length(categories),length(categories),[]);

% for c = 1:9
%   testc{c} = ones(10,1)*c;
% end
% testa = catuneven(testc,NaN);
% testa2 = reshape(testa',3,3,[])
% >>   1     4     7
%      2     5     8
%      3     6     9

if plotresults
    % Use blue to red colormap but make -1 (i.e. diagonal) black
    cmap = bwr;
    cmap(1,:) = [0 0 0];
    
    figure; imagesc(nanmean(res,3));
    colormap(cmap); caxis([-1 1]); title('mean');
    figure; imagesc(nanmedian(res,3));
    colormap(cmap); caxis([-1 1]); title('median');
end

end