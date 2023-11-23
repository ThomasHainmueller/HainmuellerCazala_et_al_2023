function [session_means, animal_means] = plot_by_animal(data,cellinfo,x,lut,varargin)
% Plot per-experiment and per-animal averages as filled, and open circles,
% respectively with consistent coloroing per animal and symbol size
% proportional to the number of cells included.
% 
% INPUT: data = column vector with raw data
%        cellinfo = struct with animal ID, session ID, etc, coming out of
%           metastatistics.
%        x = x-axis offset for this plot
%        lut = inventory of all animals with their respective proprietary
%           colors.
% 
p = inputParser;
addParameter(p,'markerScaleFct',4,@isnumeric); 
addParameter(p,'alpha',.5,@isnumeric);
addParameter(p,'LineWidth',.5,@isnumeric);

parse(p,varargin{:})

markerScaleFct = p.Results.markerScaleFct;
alpha = p.Results.alpha;
LineWidth = p.Results.LineWidth;

rng(0);

% Plot animals
animals = unique([cellinfo(:).animal]);
for a = 1:length(animals)
    idcs = find([cellinfo(:).animal] == animals(a));
    n = length(idcs);
    colID = find(lut(:,1) == cellinfo(idcs(1)).animal);
    
    %thisx = x + (rand(1,1)-.5)/2;
    thisx = x + ((a/length(animals))-0.5)/2;
    
    scat = scatter(thisx, mean(data(idcs)), n*markerScaleFct, lut(colID,2:4),'filled');
    scat.MarkerFaceAlpha = alpha;
    
    errorbar(thisx,mean(data(idcs)),std(data(idcs)),'CapSize',0,...
        'LineWidth',LineWidth,'color','k');
    
%     % SEM
%     errorbar(thisx,mean(data(idcs)),std(data(idcs))/sqrt(n),'CapSize',0,...
%         'LineWidth',1,'color','k');
    
    animal_means(a).animal = animals(a);
    animal_means(a).n = n;
    animal_means(a).data = mean(data(idcs));
end

% Plot sessions
sessions = unique([cellinfo(:).session]);
for s = 1:length(sessions)
    idcs = find([cellinfo(:).session] == sessions(s));
    n = length(idcs);
    colID = find(lut(:,1)==cellinfo(idcs(1)).animal);
    
    %thisx = x + (rand(1,1)-.5)/2;
    thisx = x + ((s/length(sessions)) - 0.5)/2;
    
    errorbar(thisx, mean(data(idcs)), std(data(idcs)), 'CapSize',0,...
        'LineWidth',LineWidth, 'Marker','o', 'MarkerSize',sqrt(n*markerScaleFct),...
        'color', [lut(colID,2:4) alpha]);
    
%     %SEM
%     errorbar(thisx, mean(data(idcs)), std(data(idcs))/sqrt(n), 'CapSize',0,...
%         'LineWidth',1, 'Marker','o', 'MarkerSize',sqrt(n*markerScaleFct), 'color', lut(colID,2:4));
    
    session_means(s).session = sessions(s);
    session_means(s).n = n;
    session_means(s).data = mean(data(idcs));
end


end