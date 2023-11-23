function data = placefields(data,method,varargin)
% Take a dataset and detect whether each cell in this dataset has a place
% field defined by the following criteria based on the dFoY traces . For
% analysis of significance, the dF values in one bin are compared to the
% distribution of dF values in the 20% bins with the lowest mean dF values.
% A place field is accepted if at least [minwidth] consecutive bins have
% significantly elevated dF values. If a place field is found that meets
% these criteria, the calculation is repeated for [nshuffle] traces in which
% dF and location data were shuffled (bootstraping) in order to give a
% p-value on how likely the place field was produced by chance. Alternative
% method adapted from Dombeck (2010,2014).
% Ftransient: Fraction of time where significant transients must be
%   present. Changed to 0.0 at 170519
if nargin<2
    method='Dombeck'; % Detection method, 'Dombeck' or 'Significance'
end

args=struct('Foutfield',0.25,'bins',0.1:0.025:2.1,'nshuffle',1000,...
    'events',true,'dFabsThresh',0.0,'dFrelThresh',7,'Ftransient',0.2,...
    'minwidth',3,'smooth',true,'minrate',0);
% 160612 Ideal Params 0.0; 7; 0.2; true based on 160501 dataset

% Overwrite default parameters if required
for pair = reshape(varargin,2,[])
    args.(pair{1}) = pair{2};
end

for c = 1:length(data.cells{1}.categories)
    ytrace = deal(data.metadata.categories{c}.y);
    ytrace = cat(2,ytrace{:}); % 170525 changed for updated DS configuration
    moving = deal(data.metadata.categories{c}.moving);
    moving = logical(cat(2,moving{:}));  
    % bootstrapped y traces.
    yrand = shuffle(ytrace(moving),args.nshuffle,50); % Default 50 frames shuffle wdw.
    for n = 1:length(data.cells)
        % Skip cells with insufficient activity
        if data.cells{n}.transientrate(c)<args.minrate
            data.cells{n}.Placefield_P(c)=1;
            data.metadata.Placefield_P(n,c)=1; %Convenience entry
            continue
        end
        
        signals=deal(data.cells{n}.categories{c}.dFoT);
        if args.events
            transients=deal(data.cells{n}.categories{c}.transientmask);
            transients=cat(2,transients{:});
            signals=cat(2,signals{:}).*transients; % 170526 adapted to new ds
        else
            signals=cat(2,signals{:});
        end
        
        if strcmp(method,'Significance')
            % Find placefield - function below.
            placefield=significant_field(signals(moving),ytrace(moving),...
                args.minwidth,args.Foutfield,args.bins,...
                args.dFabsThresh,args.Ftransient,args.smooth);
            data.cells{n}.categories{c}.placefield=placefield;
        
            % Placefields for bootstraped ytraces
            if any(placefield) % Only if a placefield has been detected.
                for ytr = args.nshuffle:-1:1 
                    bootfield(ytr)=any(significant_field(signals(moving),yrand(:,ytr),...
                args.minwidth,args.Foutfield,args.bins,...
                args.dFabsThresh,args.Ftransient,args.smooth));
                end
        
                % P-value of place field beeing genuine
                p = length(find(bootfield))/args.nshuffle;
            
                % Write to dataset:
                data.cells{n}.Placefield_P(c)=p;
            else
                % if no placefield exists.
                data.cells{n}.Placefield_P(c)=1;
            end
            data.metadata.Placefield_P(n,c)=data.cells{n}.Placefield_P(c); %Convenience entry
        end
        
        if strcmp(method,'Dombeck')
            % Find placefield - function below.
            placefield=dombeck_field(signals(moving),ytrace(moving),...
                args.minwidth,args.Foutfield,args.bins,...
                args.dFabsThresh,args.dFrelThresh,args.Ftransient,args.smooth);
            data.cells{n}.categories{c}.placefield=placefield;
        
            % Placefields for bootstraped ytraces
            if any(placefield) % Only if a placefield has been detected.
                for ytr = args.nshuffle:-1:1 
                    bootfield(ytr)=any(dombeck_field(signals(moving),yrand(:,ytr),...
                        args.minwidth,args.Foutfield,args.bins,...
                        args.dFabsThresh,args.dFrelThresh,args.Ftransient,args.smooth));
                end
        
                % P-value of place field beeing genuine
                p = length(find(bootfield))/args.nshuffle;
            
                % Write to dataset:
                data.cells{n}.Placefield_P(c)=p;
            else
                % if no placefield exists.
                data.cells{n}.Placefield_P(c)=1;
            end
            data.metadata.Placefield_P(n,c)=data.cells{n}.Placefield_P(c); %Convenience entry
        end
        %fprintf('Cell %1$i has been processed for category %2$i Place field P=%3$f\n',...
            %n,c,data.cells{n}.Placefield_P(c));
    end
end
end

function placefield = dombeck_field(dFoT,ytrace,minwidth,Foutfield,bins,...
    dFabsThresh,dFrelThresh,Ftransient,smooth)
% Function to detect place fields according to Dombeck et al., 2010,2014.
% Parameters: 'minwidth' in bins, 'dFabsThresh': minimum mean dF that must
% occur in at least one bin in the field, 'dFrelThresh' minimum 
% Distribute dF trace into bins.
placefield=false(1,length(bins)-1);
for nn = length(bins)-1:-1:1
    bindistr{nn}=dFoT((ytrace>bins(nn))&(ytrace<=bins(nn+1)));
    binmean(nn)=mean(bindistr{nn});
end
binmean(isnan(binmean))=0; % Catch cases where the bin was actually never visited

% Running mean on dF/Y trace
if smooth
    binmean=running_mean(binmean,1);
end


sortedbinmean=sort(binmean);
dFbaseline=mean(sortedbinmean(1:round(length(bins)*Foutfield)));
dFmax=max(binmean);
dFthresh=dFbaseline+0.25*(dFmax-dFbaseline); % Minimum dF for in field bins

% Identify potential place field where dF is >25% of peak dF.
for nn=length(bins)-1:-1:1
    if binmean(nn)>dFthresh
        placefield(nn)=true;
    end
end

% Eliminate all candiate fields narrower than [minwidth]
nn=1;
while nn<length(placefield)
    if nn+minwidth>length(placefield) % Treat end
        placefield(nn:end)=0;
        break
    elseif placefield(nn:nn+minwidth)
        % If true, potential PF was detected!
        PFend=nn+find(~placefield(nn:end),1,'first'); % find end of field
        if PFend>length(placefield)
            PFend=length(placefield);
        end
        % Test for absolute threshold criterion
        if any(binmean(nn:PFend)>=dFabsThresh)
            % Test for presence of significant transients F% of the time
            infielddistr=cat(2,bindistr{nn:PFend});
            % Number of nonzero entries / all entries infield
            if length(find(infielddistr>0))>Ftransient*length(infielddistr)
                nn=PFend; % Leave candidate true and jumpt to next zero value
            else
                placefield(nn:PFend)=false;
                nn=PFend; % Remove candidate and jump on to next zero value.
            end
        else
        placefield(nn:PFend)=false;
        nn=PFend; % Remove candidate and jump on to next zero value.
        end
    else
        placefield(nn)=false;
        nn=nn+1;
    end
end

% Take all detected Place fields and check whether infield values are
% [dFrelthresh] times larger than outfield.
if any(placefield) % check whether placefield exists at all
    %if ~(mean(binmean(placefield))>dFrelThresh*mean(binmean(~placefield)))
    if mean(binmean(placefield))<=dFrelThresh*mean(binmean(~placefield));
        placefield(1:end)=false;
    end
end
%plot(placefield,'r');
end

function placefield = significant_field(dFoT,ytrace,minwidth,Foutfield,bins,...
    dFabsThresh,Ftransient,smooth)
% Kernel function to placefields.m. Detect whether a place field meeting
% the criteria mentioned initially is described by a pair of dF and dY over
% time traces.
placefield=false(1,length(bins)-1);
for nn = length(bins)-1:-1:1
    bindistr{nn}=dFoT((ytrace>bins(nn))&(ytrace<=bins(nn+1)));
    binmean(nn)=mean(bindistr{nn});
end

% Running mean on dF/Y trace
if smooth
    binmean=running_mean(binmean,1);
end

sortedbinmean=sort(binmean);
dFbaseline=mean(sortedbinmean(1:round(length(bins)*Foutfield)));
dFmax=max(binmean);
dFthresh=dFbaseline+0.25*(dFmax-dFbaseline); % Minimum dF for in field bins

% Identify potential place field where dF is >25% of peak dF.
for nn=length(bins)-1:-1:1
    if binmean(nn)>dFthresh
        placefield(nn)=true;
    end
end

% Eliminate all candiate fields narrower than [minwidth]
nn=1;
while nn<length(placefield)
    if nn+minwidth>length(placefield) % Treat end
        placefield(nn:end)=0;
        break
    elseif placefield(nn:nn+minwidth)
        % If true, potential PF was detected!
        PFend=nn+find(~placefield(nn:end),1,'first'); % find end of field
        if PFend>length(placefield)
            PFend=length(placefield);
        end

        % Test for absolute threshold criterion
        if any(binmean(nn:PFend)>=dFabsThresh)
            % Test for presence of significant transients F% of the time
            infielddistr=cat(2,bindistr{nn:PFend});
            % Number of nonzero entries / all entries infield
            if length(find(infielddistr>0))>Ftransient*length(infielddistr)
                nn=PFend; % Leave candidate true and jumpt to next zero value
                % Check if infield distribution is significantly above
                % outfield.
                [~,signPF]=ranksum(cat(2,bindistr{nn:PFend})...
                    ,cat(2,bindistr{[1:nn-1,PFend+1:length(bins)-1]}));
                if ~signPF
                    placefield(nn:PFend)=false; % Remove candidate
                end
                nn=PFend; % Leave candidate true and jumpt to next zero value
            else
                placefield(nn:PFend)=false;
                nn=PFend; % Remove candidate and jump on to next zero value.
            end
        else
        placefield(nn:PFend)=false;
        nn=PFend; % Remove candidate and jump on to next zero value.
        end
    else
        placefield(nn)=false;
        nn=nn+1;
    end
end

% % Check whether value distribution in field significantly exeeds outfield.
% if any(placefield) % check whether placefield exists at all
%     [~,signPF]=ranksum(...
%         cat(2,bindistr{placefield}),cat(2,bindistr{~placefield}));
%     if ~signPF
%         placefield(1:end)=false;
%     end
% end
end
