function quotients = famnov_meta(dataset,categories,varargin)
% Plot activity differences between novel and familiar. Designed for
% interneuron data. Plot quotient of fam and nov activity separately for
% each run and speedmaps of activity.

if nargin < 2
    categories = [1,2];
end

args=struct('nruns',0,'xlim',[1 15],'ylim',[0 .5]); 
    % nruns: specify, if the number of runs is inhomogenous between
    %   dataset. By default, nruns of the first dataset is used.

% Overwrite default parameters if required
for pair = reshape(varargin,2,[])
    %args.(pair{1}) = pair{2};
    if isfield(args,pair{1})
        args.(pair{1})=pair{2};
    else
        error('Input argument does not exist\n')
    end
end

quotients = []; % response indices
dFodY = [];

% Get Stimulus PSTH for all cells throughout datasets.
for d=1:numel(dataset)
    ds=dataset{d};
    
    if ~args.nruns
        args.nruns = length(ds.cells{1}.categories{categories(1)}.dFoT);
    end
    
    for n=1:length(ds.cells)
        for r=1:args.nruns
            thisquotient(r) = nanmean(ds.cells{n}.categories{categories(2)}.dFoT{r})/...
                nanmean(ds.cells{n}.categories{categories(1)}.dFoT{r});
            thisdFodY(:,r,1) = ds.cells{n}.categories{categories(1)}.dFodY{r};
            thisdFodY(:,r,2) = ds.cells{n}.categories{categories(2)}.dFodY{r};
        end
        quotients(:,end+1) = thisquotient;
        dFodY(:,:,:,end+1) = thisdFodY;
    end
end

dFodY(:,:,:,1) = [];

% Display part
figure;
for r = 1:args.nruns
    subplot(1,args.nruns,r)
        hold on
        xlim(args.xlim);
        ylim(args.ylim);
        
        % Familiar (black - blue)
        plot(nanmean(dFodY(:,r,1,:),4),'k');
        plot(nanmean(dFodY(:,r,1,:),4)+...
            nanstd(dFodY(:,r,1,:),0,4)/sqrt(size(dFodY,4)),'b');    
        plot(nanmean(dFodY(:,r,1,:),4)-...
            nanstd(dFodY(:,r,1,:),0,4)/sqrt(size(dFodY,4)),'b');    
        
        % Novel (red - magenta)
        plot(nanmean(dFodY(:,r,2,:),4),'r');
        plot(nanmean(dFodY(:,r,2,:),4)+...
            nanstd(dFodY(:,r,2,:),0,4)/sqrt(size(dFodY,4)),'m');    
        plot(nanmean(dFodY(:,r,2,:),4)-...
            nanstd(dFodY(:,r,2,:),0,4)/sqrt(size(dFodY,4)),'m');      
end
end
