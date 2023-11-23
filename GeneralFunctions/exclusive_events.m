function data = exclusive_events(data,varargin)
% Take a dataset obtained with get_lineartracks and corrected category
% settings and apply the standard working procedures.
args=struct('display',1);

% Overwrite default parameters if required
for pair = reshape(varargin,2,[])
    args.(pair{1}) = pair{2};
end

for c=1:length(data.metadata.categories)
    
    nsomatic = 0;
    ndendritic = 0;
    nsomatodendritic = 0;
    
    for r = 1:length(data.metadata.categories{c}.y)
        % Collect all data for this run from soma and dendrites
        somamask = data.soma{1}.categories{c}.transientmask{r};
        for d = 1:length(data.dendrites)
            denmask(d,:) = data.dendrites{d}.categories{c}.transientmask{r};
        end
        alldenmask = max(denmask,[],1);        
        
        % Prepare arrays for this run
        somaexclusive = false(1,length(somamask));
        denexclusive = false(1,length(somamask));
        somatodendritic = false(1,length(somamask));
        
        % 1st round: classify traces with somatic involvement
        somatimes = [find(diff(somamask)==1)+1,length(somamask)];
        for t = 1:length(somatimes)-1
            period = somatimes(t):somatimes(t+1)-1;
            % Does this somatic transient overlap with any dendritic one?
            if any(and(somamask(period),alldenmask(period)))
                somatodendritic(period) = somamask(period);
            else
                somaexclusive(period) = somamask(period);
            end
        end
        
        % 2nd round: get dendrite-exclusive events
        for d = 1:length(data.dendrites)
            thisdentimes = [find(diff(denmask(d,:))==1)+1,length(somamask)];
            for t = 1:length(thisdentimes)-1
                period = thisdentimes(t):thisdentimes(t+1)-1;
                if any(and(denmask(d,period),...
                        max(cat(1,somatodendritic(period),somaexclusive(period),denexclusive(period)),[],1)))
                    continue % Transient overlaps with a detected one
                else
                    denexclusive(period) = denmask(d,period);
                end
            end
        end
        
        data.metadata.categories{c}.soma_events{r} = somaexclusive;
        data.metadata.categories{c}.dendrite_events{r} = denexclusive;
        data.metadata.categories{c}.somatodendritic_events{r} = somatodendritic;
        
        nsomatic = nsomatic + length(find(diff(somaexclusive)==1));
        ndendritic = ndendritic + length(find(diff(denexclusive)==1));
        nsomatodendritic = nsomatodendritic + length(find(diff(somatodendritic)==1));
        
        if args.display
            figure; hold on;
            
            % Plot somatic transients
            somatrace = data.soma{1}.categories{c}.dFoT{r};
            somatrans = somatrace.*data.soma{1}.categories{c}.transientmask{r};
            somatrans(somatrans==0) = NaN;
            plot(somatrace,'b');
            plot(somatrans,'r');
            
            % Plot dendritic transients
            for d=1:length(data.dendrites)
                dentrace = data.dendrites{d}.categories{c}.dFoT{r};
                dentrans = dentrace.*data.dendrites{d}.categories{c}.transientmask{r};
                dentrans(dentrans==0) = NaN;
                plot(dentrace-d,'k');
                plot(dentrans-d,'r');
            end
            
            % Plot classification results
            plot(somaexclusive*(d+2)-(d+1),'b');
            plot(denexclusive*(d+2)-(d+1),'k');
            plot(somatodendritic*(d+2)-(d+1),'c');
        end
        clear somamask denmask
    end
    
    data.metadata.nsomatic(c) = nsomatic;
    data.metadata.ndendritic(c) = ndendritic;
    data.metadata.nsomatodendritic(c) = nsomatodendritic;
    
end

end