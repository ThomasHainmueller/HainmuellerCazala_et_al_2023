function [pPost,pPostMax,SqErr] = BayesWrapper(data, varargin)
% Get traces from the selected categories of a dataset and perform bayesian
% decoding of the trajectory using the traces from specified cells
% Optional arguments:
% categories:	int, default: use all categories of the dataset
% cells:        int, specifiy the cells used for en- and decoding
% nsegments:	int, decoding is performed on segments omited during
%           	encoding. Specify the number of segements.
% moving_only:  bool, default true.
% transients_only:   bool, whether to use only detected transients or the
%                   entire fluorescence trace.
% spbins:         double, spatial bins for discretizing the spatial maps
%
% Note to self: Positional- and Contextual information are the respective
% marginals of the two decoded probability distributions.
%% Input handling
p = inputParser;
addParameter(p,'categories',[],@isnumeric);
addParameter(p,'cells',[],@isnumeric);
addParameter(p,'nsegments',10,@isnumeric);
addParameter(p,'moving_only',true,@islogical);
addParameter(p,'transients_only',false,@islogical);
addParameter(p,'spbins',0.1:0.025:2.1,@isnumeric);

parse(p,varargin{:})

nsegments = p.Results.nsegments;
moving_only = p.Results.moving_only;
transients_only = p.Results.transients_only;
spbins = p.Results.spbins;

if ~isempty(p.Results.categories)
    categories = p.Results.categories;
else
    categories = 1:length(data.metadata.categories);
end

if ~isempty(p.Results.cells)
    cells = p.Results.cells;
else
    cells = 1:length(data.cells);
end

%% Get the animal's position, moving intervals and traces
for c = categories
    pos{c} = deal(data.metadata.categories{c}.y(:));
    pos{c} = cat(2,pos{c}{:});
    
    % Make bins for interleaved decoding.
    segments{c} = round(linspace(1,size(pos{c},2),nsegments+1));
    
    if moving_only
        moving{c} = deal(data.metadata.categories{c}.moving(:));
        moving{c} = cat(2,moving{c}{:});
    else
        moving{c} = true(1,size(pos{c},2));
    end
    
    % Get traces
    for n = flip(cells)
        temptr = deal(data.cells{n}.categories{c}.dFoT(:));
        temptr = cat(2,temptr{:});
        
        % When it sees negative values, the bayesian decoder becomes very angry
        % and throws complex numbers at you, hence:
        temptr(temptr<0) = 0;
        
        if transients_only
            mask = deal(data.cells{n}.categories{c}.transientmask{:});
            mask = cat(2,mask{:});
            temptr = temptr.*mask;
        end
        
        ctraces(:,n) = temptr;
    end
    traces{c} = ctraces;
    clear ctraces
end

%% Encoding and decoding
for b = 1:nsegments
    for c = categories
        iencode = true(1,size(traces{c},1));
        iencode(1,segments{c}(b):segments{c}(b+1)) = false; % Remove decoding bin
        %idecode{c} = ~iencode;
        
        iencode = and(iencode,moving{c});
        %idecode{c} = and(idecode{c},moving{c});
        
        % Make placemaps aka encoding
        for n = flip(cells)
            ratemaps(n,:,c) = SBdiscretize(traces{c}(iencode,n),...
            	pos{c}(iencode),spbins);

            % Use Gaussian filter to smooth the space plot.
%             nanvals = isnan(rawdFoY);
%             rawdFoY(nanvals)=0;
%             rawdFoY = conv(rawdFoY,filter,'same');
%             rawdFoY(nanvals)=NaN;
%             thisdFoY = cat(2,thisdFoY,rawdFoY);
        end
    end
    % Linearize ratemaps
    ratemaps = reshape(ratemaps,size(ratemaps,1),size(ratemaps,2)*length(categories));
    ratemaps = ratemaps';
    
    % Decoding
    for c = categories
        idecode = false(1,size(traces{c},1));
        idecode(1,segments{c}(b):segments{c}(b+1)) = true; % Keep Decoding bin
        idecode = and(idecode,moving{c});
        
%         pPost{c}(idecode,1:size(ratemaps,1)) = bayes_decode(ratemaps,...
%             traces{c}(idecode,:)',1/data.metadata.categories{c}.acquisition_rate(1));
        [~,pPost{c}(idecode,1:size(ratemaps,1))] = PVcorr_decode(ratemaps,...
            traces{c}(idecode,:)');
    end
    clear ratemaps
end

%% Compute mean square error of decoded position.
end
