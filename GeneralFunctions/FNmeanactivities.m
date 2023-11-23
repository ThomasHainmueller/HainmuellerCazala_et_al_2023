function data = FNmeanactivities(data, category)
% For 'PVcells' type datasets. Preliminary function to get the mean DF/F
% activity for familiar and novel context for each cell.
if nargin<2
    category = 1;
end

for n=1:length(data.cells)
    for r = 1:length(data.metadata.categories{category}.x)
        cattrace = data.metadata.categories{category}.x{r};
        cattrace = medfilt1(cattrace, 50); %Median filtering over 50 bins to xclude artefacts
        thisfam = cattrace<1; % Works with standard fam =.5 and nov=2.0 V in x-coordinate
        thisnov = cattrace>=1;
        
        thistrace = data.cells{n}.categories{category}.dFoT{r};
        
        data.cells{n}.categories{category}.FAMmean(r) = nanmean(thistrace(thisfam))+1;
        data.cells{n}.categories{category}.NOVmean(r) = nanmean(thistrace(thisnov))+1;
    end
end
    
end