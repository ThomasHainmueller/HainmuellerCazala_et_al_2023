function data = PlacefieldsPosNeg(data)
% Response to reviewer comment to Hainmueller, Cazala, et al., 
%
% Thomas Hainmueller, Buzsakilab, 2023

bins = .1:.025:1.9;

%% Extract position and calcium traces, interpolate to 10 Hz framerate (100 ms bins)
for c = 1%:length(data.metadata.categories)
    pos = cat(2,data.metadata.categories{c}.y{:});
    mov = cat(2,data.metadata.categories{c}.moving{:}) > 0;
    
    % Remove immobility periods
    pos = pos(mov);
    
    for n = 1:length(data.cells)
        trace = cat(2,data.cells{n}.categories{c}.dFoT{:});
        trace = trace(mov);
        
        placefield = find_fields(trace, pos, bins);
        
        %figure; hold on; plot(placefield); plot(SBdiscretize(trace, pos, bins, 'median'));
        
        data.cells{n}.categories{c}.placefield = placefield;
        data.cells{n}.nPosFields(c) = length(find([0 diff(placefield)] == 1 & placefield == 1));
        data.cells{n}.nNegFields(c) = length(find([0 diff(placefield)] == -1 & placefield == -1));
    end
end

end

function placefields = find_fields(trace, pos, bins)
% Find 'placefields' as continuous series of locations where mean activity
% distribution deviates significantly above (positive) or below (negative)
% distribution across all bins.

minwidth = 4; % In bins
sigthr = 0.05; % significance threshold
minoccup = 10; % Min number of samples per bin

binhist = SBdiscretize(trace,pos,bins,'distribution');

% Create groups for Kruskal-Wallis (below)
posIDs = discretize(pos,bins);
posIDs = SBdiscretize(posIDs,pos,bins,'distribution');

% Perform Kruskal-wallis test for uneven distribution of fluorescence
[pKW,~,stats] = kruskalwallis(cat(2,binhist{:}),cat(2,posIDs{:}),'off');
%multcompare(stats);

if pKW >= sigthr
    placefields = zeros(1,length(binhist));
    return
end

% Adjust significance threshold (bonferoni correction)
sigthr = sigthr/length(binhist);

% Find bins that are significantly above (1) or below (-1) the rest
for b = length(binhist):-1:1
    if length(find(~isnan(binhist{b}))) < minoccup
        placefields(b) = 0;
        continue
    end
        
    otherhist = cat(2,binhist{setdiff(1:length(binhist),b)});
    
    [~,p] = ranksum(otherhist,binhist{b});
    ispos = median(binhist{b}) > median(otherhist);
    
    if p < sigthr & ispos
        placefields(b) = 1;
    elseif p < sigthr & ~ispos
        placefields(b) = -1;
    else
        placefields(b) = 0;
    end
end

% Remove 'fields' that are shorter than 'minwidth'
% Positive fields
F = find(diff([0,placefields>0,0]));
[~,X] = find(diff(F)<minwidth);
for k = X
    placefields(F(k):F(k+1)-1) = 0; 
end

% Negative fields
F = find(diff([0,placefields<0,0]));
[~,X] = find(diff(F)<minwidth);
for k = X
    placefields(F(k):F(k+1)-1) = 0; 
end


% b = 1;
% 
% while b < length(binhist)
%     % Treat ends
%     if length(binhist)-b <= minwidth
%         if ~placefields(b:end) == placefields(b)
%             placefields(b:end) = 0;
%         end
%         b = length(binhist);
%         
%     % Skip if no placefield in this bin
%     elseif placefields(b) == 0
%         b = b+1;
%     
%     % Check if placefield meets length criteria
%     elseif all(placefields(b:b+minwidth-1) == placefields(b))
%         nextout = find(placefields(b:end) ~= placefields(b),1,'first');
%         if ~isempty(nextout)
%             b = b + nextout -1;
%         else
%             b = length(binhist);
%         end
%     else
%         placefields(b) = 0;
%         b = b+1;
%     end
% end

end

