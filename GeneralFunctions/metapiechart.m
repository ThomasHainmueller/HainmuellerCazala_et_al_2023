function results = metapiechart(dataset, parameter, varargin)
% Takes metadataset, finds for TWO categories the number of cells that
% fulfill a certain parameter in either or both of the categories. Returns
% a 4 x 'number of ds in metads' table with column order A, B, both, none
% and plots a pie chart of the mean results.
% paramters: 'active', 'PF', 'SI'

args=struct('active',0,'SI',1.0,'PF',1.0,'categories',...
    [1,2], 'combination','or'); 
    % SI/PF significance levels of spatial info/placefields from bootstrap

% Overwrite default parameters if required
for pair = reshape(varargin,2,[])
    %args.(pair{1}) = pair{2};
    if isfield(args,pair{1})
        args.(pair{1})=pair{2};
    else
        error('Input argument does not exist\n')
    end
end

for ds = length(dataset):-1:1
    results(3,ds) = length(findcells(dataset{ds}, args.categories, args.active,...
        args.SI, args.PF ,'and'));
    results(1,ds) = length(findcells(dataset{ds}, args.categories(1), args.active,...
        args.SI, args.PF ,'or')) - results(3,ds);
    results(2,ds) = length(findcells(dataset{ds}, args.categories(2), args.active,...
        args.SI, args.PF ,'or')) - results(3,ds);
    
    switch parameter
        case 'active'
            results(4,ds) = length(findcells(dataset{ds}, args.categories, 0,...
                args.SI, args.PF ,'or')) - sum(results(1:3,ds));
        case 'PF'
            results(4,ds) = length(findcells(dataset{ds}, args.categories,...
                args.active, args.SI, 1 ,'or')) - sum(results(1:3,ds));
        case 'SI'
            results(4,ds) = length(findcells(dataset{ds}, args.categories,...
                args.active, 1, args.SI ,'or')) - sum(results(1:3,ds));
    end
end

%figure; pie(sum(results,2), {'A','B','A and B','none'});
figure; pie(sum(results,2));

results = results';
end