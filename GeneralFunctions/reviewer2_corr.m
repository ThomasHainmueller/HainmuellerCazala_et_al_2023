function data = reviewer2_corr(data,sigma)
% Computes the trial-to-trial correlation within each category and appends
% the category means to the 'cell' data.
% Also computes the 1st vs. 2nd half correlation of the runs and appends
% them in the same manner.

if nargin<2
    sigma=1; % Gaussian function with sigma = 1 frame, set 0 to disable
end

if sigma
    filter = fspecial('gaussian',[6*sigma,1],sigma);
end

% Find the correlations between categories for each cell individually.
for n=length(data.cells):-1:1
    for c=length(data.cells{n}.categories):-1:1
        % Get dFoY traces for this cell and category
        thisdFoY=deal(data.cells{n}.categories{c}.dFoY);
        thisdFoY=cat(2,thisdFoY{:});

        % Calculate trial-to-trial variability. Smoothing of traces may not be
        % advisable here
        trialvar = corrcoef(thisdFoY,'rows','pairwise');
        % set autocorrelations to NaN
        trialvar(1:(size(trialvar,1)+1):size(trialvar,1)*size(trialvar,2)) = NaN;
        data.cells{n}.trialvar(c)=...
            nanmean(reshape(trialvar,size(trialvar,1)*size(trialvar,2),1));
        
        % Calculate 1st vs. 2nd half correlations within each category
        if sigma
            nanvals = isnan(thisdFoY);
            thisdFoY(nanvals)=0;
            for r=1:size(thisdFoY,2)
                thisdFoY(:,r) = conv(thisdFoY(:,r),filter,'same');
            end
            thisdFoY(nanvals)=NaN;
        end
        
        midpoint = round(size(thisdFoY,2)/2); % To split the trials in early and late half
        stability = corrcoef(nanmean(thisdFoY(:,1:midpoint),2),...
            nanmean(thisdFoY(:,midpoint+1:end),2),'rows','complete');
        data.cells{n}.sessionstab(c)= stability(1,2);
    end
end
end