function [indices,mask] = findcells(data,categories,active,pSI,pPF,combination,region)
% Get indices of cells in a dataset that pass the set threshold for
% activity and significance of spatial information.
%
% CAUTION: BEHAVIOR OF THIS FUNCTION WAS CHANGED 221204. Before, the 'or'
% option returned cells that met each of the criteria in at least one
% condition (e.g. activity on day one and SI on day 2). Going forward, it
% will only return cells that meet ALL of the criteria in at least one of
% the conditions (e.g. activity and SI threshold in the same category.
% The 'and' option always returned cells that meet ALL criteria in ALL
% selected categoires!

if nargin<7
    region = '';
end
if nargin<6
    combination='and';
end
if nargin<5
    pPF=1; % Significance level for having a placefield
end
if nargin<4
    pSI=0.05; % default singificance of spatial information>chance
end
if nargin<3
    active = 1/60; % default: >1 transient/min
end

% activecells=false(1,length(data.cells));
% SIcells=false(1,length(data.cells));
% PFcells=false(1,length(data.cells));
inregioncells = false(length(data.cells),1);

% ADDED 221210, exclude categories that are not part of the particular DS
selcategories = intersect(categories,1:length(data.metadata.categories));

if isempty(selcategories)
    mask = false(1,length(data.cells));
    
elseif length(selcategories) ~= length(categories) && strcmp(combination,'and')
    % Return empty if selection is 'and' and not all categories exist
    mask = false(1,length(data.cells));
    
else
    for n = length(data.cells):-1:1
        
        % Find active cells
        transientrates(n,:)=data.cells{n}.transientrate(selcategories);
        SIpvals(n,:) = data.cells{n}.spatial_P(selcategories);
        PFpvals(n,:) = data.cells{n}.Placefield_P(selcategories);
        %transientrates=cat(2,transientrates{:});
        %transientrates=transientrates(categories);
        %     if strcmp(combination,'and') && all(transientrates>=active)
        %         activecells(n)=true;
        %     elseif strcmp(combination,'or') && any(transientrates>=active)
        %         activecells(n)=true;
        %     end
        
        % Find spatially modulated cells
        %     SI=data.cells{n}.spatial_P(categories);
        %     %SI=cat(2,SI{:});
        %     %SI=SI(categories);
        %     if strcmp(combination,'and') && all(SI<=pSI)
        %         SIcells(n)=true;
        %     elseif strcmp(combination,'or') && any(SI<=pSI)
        %         SIcells(n)=true;
        %     end
        
        % Find cells that have a defined placefield
        %     PF=data.cells{n}.Placefield_P(categories);
        %     %PF=cat(2,PF{:});
        %     %PF=PF(categories);
        %     if strcmp(combination,'and') && all(PF<=pPF)
        %         PFcells(n)=true;
        %     elseif strcmp(combination,'or') && any(PF<=pPF)
        %         PFcells(n)=true;
        %     end
        
        % Select cells that lie in a given region (e.g. superficial)
        if region
            if isfield(data.cells{n},'region')
                if strcmp(data.cells{n}.region, region)
                    inregioncells(n,1)=true;
                end
            end
        else
            inregioncells(n,1)=true;
        end
    end
    
    selection(:,:,1) = transientrates >= active;
    selection(:,:,2) = SIpvals <= pSI;
    selection(:,:,3) = PFpvals <= pPF;
    
    if strcmp(combination,'and')
        meetcriteria = all(selection,[2 3]);
    elseif strcmp(combination,'or')
        meetcriteria = any(all(selection,3),2);
    end
    
    mask = meetcriteria & inregioncells;
end

% Rest of code uses row-vectors unfortunately
mask = reshape(mask,1,[]);
indices=find(mask);
end