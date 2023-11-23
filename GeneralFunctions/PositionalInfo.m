function data = PositionalInfo(data)
% Get positional information score (Olypher and Fenton, 2003; Wilent, Nitz 2007) 
% for imaging dataset.
%
% Thomas Hainmueller, Buzsakilab, 2023

bins = .1:.025:1.9;
nshuffles = 1000;
rng(0);

%% Extract position and calcium traces, interpolate to 10 Hz framerate (100 ms bins)
for c = 1%:length(data.metadata.categories)
    for r = 1:length(data.metadata.categories{c}.y)
        pos{r} = data.metadata.categories{c}.y{r};
        mov = data.metadata.categories{c}.moving{r} > 0;
        
        % Remove immobility periods
        pos{r} = pos{r}(mov);

        % Generate timestamps for original framerate (do not use ft, as
        % they also contain immobility periods)
        origts = (0:length(pos{r})-1)/data.metadata.categories{c}.acquisition_rate(r);
        pos{r} = interp1(origts,pos{r},0:.1:origts(end));
        
        for n = 1:length(data.cells)
            thistr = data.cells{n}.categories{c}.dFoT{r};
            thistr = thistr(mov);
            
            % Remove <0 values from calcium traces
            thistr(thistr<0)=0;
            tr(n,:) = interp1(origts,thistr,0:.1:origts(end));
        end
        trcs{r} = tr;
        clear tr
    end
    
    position = cat(2,pos{:});
    traces = cat(2,trcs{:});
    
    % Constrain position to bin limits
    position(position<bins(1)) = bins(1);
    position(position>bins(end)) = bins(end);
    
    % Compute positional information
    for n = 1:length(data.cells)
        % Replace extreme outliers in activity distribution with maximum
        % non-outlier value (outliers can only be positive, <0 removed).
        trace = traces(n,:);
        [~,isoutlier] = tukeyOutlierRemoval(trace',3);
        trace(isoutlier) = max(trace(~isoutlier));
        
        data.cells{n}.positional_info(c) = getMutualInfo(position,trace,bins);
        
        % Bootstrap here
        randinfo = NaN(1,nshuffles);
        parfor sh = 1:nshuffles
            randinfo(sh) = getMutualInfo(...
                position,circshift(trace,round(rand*length(trace))),bins);
        end
        
        data.cells{n}.positional_info_rand{c} = randinfo;
        data.cells{n}.positional_info_r(c) = data.cells{n}.positional_info_rand{c}(1);
        data.cells{n}.positional_P(c) = length(find(...
            randinfo >= data.cells{n}.positional_info(c)))/nshuffles;
    end
end
end


function Ipos = getMutualInfo(position,trace,bins)
% Compute the mean positional info value (i.e. mutual information between
% position and activity)
% 
% Ipos(xi) = P(activitylevel|xi) * log2(P(activitylevel|xi)/P(activitylevel))


% Discretize activity levels; Somewhat backwards but need this to create
% the bins
[~,trbins] = discretize(trace,10);
histAP = hist3(cat(1,trace,position)','Edges',{trbins bins});
histAP = histAP(1:end-1,1:end-1); % Remove extra bins added by function (empty).

% Compute Positional info for each bin
% Ipos(xi) = P(activitylevel(k)|xi) * log2(P(activitylevel(k)|xi)/P(activitylevel))
for xi = size(histAP,2):-1:1
    Pkxi = histAP(:,xi)./sum(histAP(:,xi));
    Pk = sum(histAP,2)./sum(histAP,'all');
    Ipos(xi) = nansum( Pkxi .* log2(Pkxi./Pk) );
end

Ppos = sum(histAP,1)/sum(histAP,'all');
Ipos = sum(Ipos.*Ppos);

end