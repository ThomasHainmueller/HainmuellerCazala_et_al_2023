function resmatrix = metacorr(dataset, varargin)
% Rather specific. Takes metads and creates the correlation matrix in the
% order specified. Each row is for the cells that fullfill specified
% criteria for the respective category.

args=struct('active',0,'SI',1.0,'PF',1.0,'region','','categories',...
    1:length(dataset{1}.metadata.categories), 'combination','or',...
    'parameter','Pearson_r'); 
    % SI/PF significance levels of spatial info/placefields from bootstrap
    % Categories can be used to specify the order of the resulttable
    % parameter: Works with Random_Pearson_r as well (for comparisons)

% Overwrite default parameters if required
for pair = reshape(varargin,2,[])
    %args.(pair{1}) = pair{2};
    if isfield(args,pair{1})
        args.(pair{1})=pair{2};
    else
        error('Input argument does not exist\n')
    end
end

n=1;

for c = args.categories
    thismatrix = metastatistics(dataset,'cell',args.parameter,'categories',c,...
        'active',args.active,'PF',args.PF,'SI',args.SI,'region',args.region,...
        'rescategories',1:length(args.categories));
    close();
    
    results{n} = thismatrix(:,c,args.categories); %vector
    %resmatrix(n,:) = squeeze(nanmean(results{c},1));
    %medmatrix(n,:) = squeeze(nanmedian(results{c},1));
    ncells(c) = size(thismatrix,1);
    n=n+1;
end

% reorganize the results into a matrix padded with NaNs
resmatrix = NaN(max(ncells),c,c);

for c = 1:length(args.categories)
    resmatrix(1:size(results{c},1),c,1:length(args.categories)) = results{c};
end

figure; imagesc(squeeze(nanmean(resmatrix,1)));
    colormap('jet'); caxis([.1 .8]); title('mean');
figure; imagesc(squeeze(nanmedian(resmatrix,1))); 
    colormap('jet'); caxis([.1 .8]); title('median');
%figure; imagesc(resmatrix); colormap('jet'); caxis([.1 .8]); title('mean');
%figure; imagesc(medmatrix); colormap('jet'); caxis([.1 .8]); title('median');
end