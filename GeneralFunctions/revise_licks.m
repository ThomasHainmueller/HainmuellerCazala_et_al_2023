function [dataset] = revise_licks(dataset, varargin)
% Rerun the lick detection. Should be run from the folder in which the
% '*.ephys' files are stored.

args=struct('folder',[],'display',false,'bins',.1:.05:2.1);

% Overwrite default parameters if required
for pair = reshape(varargin,2,[])
    if isfield(args,pair{1})
        args.(pair{1})=pair{2};
    else
        error('Input argument does not exist\n')
    end
end

for c = 1:length(dataset.metadata.categories)
    for run = 1:length(dataset.metadata.categories{c}.licks)
        fn = dataset.metadata.categories{c}.filename{run};
        
        if args.display
            ephys = get_behaviourdata(fn,3,1000,true); % for debugging only
        else
            ephys = get_behaviourdata(fn);
        end
        
        dataset.metadata.categories{c}.licks{run} = ephys(:,3);
    end
end