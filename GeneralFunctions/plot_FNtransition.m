function plot_FNtransition(data,n,signal,SDscale)
% Plot calcium traces aligned to the familiar to novel transition for
% individual cells or boutons.
% 
% Thomas Hainmueller, 2023
if nargin<4
    SDscale = 10;
end
if nargin<3
    signal = 'dFoT';
end

colors = {[1 .5 .5],[1 .7 .7]};

fr = data.metadata.categories{1}.acquisition_rate{1};
nruns = min(length(data.metadata.categories{1}.y),...
    length(data.metadata.categories{2}.y));
spacing = cat(2,data.cells{n}.categories{1}.(signal){:},...
    data.cells{n}.categories{1}.(signal){:});
spacing = SDscale*nanstd(spacing);

figure; hold on;

for r = 1:nruns
    Ftr = data.cells{n}.categories{1}.(signal){r};
    Ntr = data.cells{n}.categories{2}.(signal){r};
    
    Fts = (-length(Ftr):-1)/fr;
    Nts = (0:length(Ntr)-1)/fr;
    
    plot(Fts,Ftr - r*spacing,'color',colors{1});
    plot(Nts,Ntr - r*spacing,'color',colors{2});
end
    
ylabel('DF/F');
xlabel('Time (s)');

end