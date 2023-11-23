function [ds] = remove_duplicates(ds,corrthr,varargin)
% Remove putative duplicate elements (cells or boutons) from dataset based
% on the element-by-element correlation structure. Elements which have
% calcium signals exceeding a threshold of 'corrthr' are excluded from the
% dataset.

p = inputParser;
addParameter(p,'actThr',0,@isnumeric); % Activitiy threshold
addParameter(p,'categories',1,@isnumeric); % Plot output to command line, needs to be ROW vector (i.e. [1 2])
addParameter(p,'smoothFct',1,@isnumeric); % window width, 1 = disabled, 5 = matlab standard
addParameter(p,'signal','dFoT',@ischar); % signal to extract, 'dFoT','dZoT','dFoT_fneu'

parse(p,varargin{:})

actThr = p.Results.actThr;
categories = p.Results.categories;
smoothFct = p.Results.smoothFct;
signal = p.Results.signal;

% Get the concatenated dFoT traces of all cells
for n = length(ds.cells):-1:1
    thistr = [];
    for c = categories
        thistr = cat(2,thistr,ds.cells{n}.categories{c}.(signal){:}); % 1st category only
    end
    traces(n,:) = smooth(thistr,smoothFct); %cat(2,thistr{:});
    
    if actThr > 0
        inactive(n) = ds.cells{n}.transientrate(1)<actThr;
    else
        inactive(n) = false;
    end
end

corrmtrx = corrcoef(traces');
corrmtrx(find(eye(size(corrmtrx)))) = NaN;
%figure; imagesc(corrmtrx)

[maxcorr,mostsimilar] = max(corrmtrx);

% Iteratively remove exceedingly correlated elements
isduplicate = false(1,length(ds.cells));
% maxcorr_orig = maxcorr; % Keep unaltered copy

while max(maxcorr) > corrthr
    [~,idx] = max(maxcorr);
    isduplicate(idx) = true;
    corrmtrx(:,idx) = NaN;
    corrmtrx(idx,:) = NaN;
    maxcorr = max(corrmtrx);
end

%maxcorr_orig(inactive) = NaN;
isduplicate(inactive) = false;

% Remove duplicate cells
ds.cells = ds.cells(~isduplicate);

% Add similarity values to cells
% for n = 1:length(mds{nds}.cells)
%     ds.cells{n}.maxcorr = maxcorr_orig(n);
%     ds.cells{n}.mostsimilar = mostsimilar(n);
%     ds.cells{n}.isduplicate = isduplicate(n);
% end


end