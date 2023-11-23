function [indices,mask] = findcells_meta(mds,categories,active,pSI,pPF,combination,region)
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

for ds = 1:length(mds)
    [indices{ds} mask{ds}] = findcells(mds{ds},categories,...
        active,pSI,pPF,combination,region);
end

mask = cat(2,mask{:});
end