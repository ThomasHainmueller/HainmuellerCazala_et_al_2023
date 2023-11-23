function results = analyze_dataset_meta(metads, parameter, varargin)
% Take a metadataset and extract values for a selected parameter.
%
% VARIABLES:
% metads = cell with dataset structures
% parameter = string, must be present in all datasets
% compunit = string, ['cell','dataset'], unit of comparison, can be all
%   individual cells or the mean of the dataset for this parameter
% section = string, ['cells', 'metadata']
% selection = bool, use only data from selected cells as specified. See
%   findcells.m for parameter specification.
% categories = int, specify for metads with unequal category counts

args=struct('compunit','cell','section','cells','selection',0,'categories',...
    1:length(metads{1}.metadata.categories),'active',0,'pSI',1,'pPF',1,...
    'combination','or');

% Overwrite default parameters if required
for pair = reshape(varargin,2,[])
    if isfield(args,pair{1})
        args.(pair{1})=pair{2};
    else
        error('Input argument does not exist\n')
    end
end

% Data extraction
results = [];
% TODO: Implement cell selection acc to criteria
for nds = 1:length(metads)
	if strcmp(args.section,'cells')
        subselection = findcells(metads{nds}, args.categories,...
            args.active, args.pSI, args.pPF, args.combination);
    else
        subselection = 1; % For non
    end
    
    for n = 1:length(metads{nds}.(args.section)(subselection))
        results(end+1,:)=...
            metads{nds}.(args.section){n}.(parameter)(args.categories); % Very academic, reads e.g. to "meta{1}.cells{20}.transientrate(1:4)"
    end
end
end