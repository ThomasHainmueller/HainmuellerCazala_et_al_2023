function [session_means, animal_means] = scatter_by_animal(data,cellinfo,lut,varargin)
% Plot per-experiment and per-animal averages as filled, and open circles,
% respectively with consistent coloroing per animal and symbol size
% proportional to the number of cells included. This plot is for paired
% values (e.g. bootstrap vs. actual or fam vs. nov) and comes with vertical
% and horizontal errorbars.
% 
% INPUT: data = Nx2 column vector with raw data
%        cellinfo = struct with animal ID, session ID, etc, coming out of
%           metastatistics.
%        lut = inventory of all animals with their respective proprietary
%           colors.
% 
p = inputParser;
addParameter(p,'markerScaleFct',2,@isnumeric); 
addParameter(p,'alpha',.5,@isnumeric);
addParameter(p,'LineWidth',.5,@isnumeric);

parse(p,varargin{:})

markerScaleFct = p.Results.markerScaleFct;
alpha = p.Results.alpha;
LineWidth = p.Results.LineWidth;

rng(0);

% Plot sessions
sessions = unique([cellinfo(:).session]);
for s = 1:length(sessions)
    idcs = find([cellinfo(:).session] == sessions(s));
    %n = length(idcs);
    n = length(find(~isnan(data(1,idcs)) | ~isnan(data(2,idcs)))); % Cave, this assumes unpaired data - check on input side!
    colID = find(lut(:,1)==cellinfo(idcs(1)).animal);
    
    errorbar(nanmean(data(1,idcs)), nanmean(data(2,idcs)),...
        nanstd(data(2,idcs)),nanstd(data(2,idcs)),...
        nanstd(data(1,idcs)),nanstd(data(1,idcs)),... 
        'CapSize',0,'LineWidth',LineWidth, 'Marker','o',...
        'MarkerSize',sqrt(.1+n*markerScaleFct), 'color', [lut(colID,2:4) alpha]);  % +.1 if n == 0
    
%     % SEM
%     errorbar(mean(data(2,idcs)), mean(data(1,idcs)),...
%         std(data(1,idcs))/sqrt(n),std(data(1,idcs))/sqrt(n),...
%         std(data(2,idcs))/sqrt(n),std(data(2,idcs))/sqrt(n),...
%         'CapSize',0,'LineWidth',1, 'Marker','o',...
%         'MarkerSize',sqrt(n*markerScaleFct), 'color', lut(colID,2:4));
    
    session_means(s).session = sessions(s);
    session_means(s).n = n;
    session_means(s).data = nanmean(data(:,idcs),2);
    session_means(s).SEM = nanstd(data(:,idcs),[],2);
end

% Plot animals
animals = unique([cellinfo(:).animal]);
for a = 1:length(animals)
    idcs = find([cellinfo(:).animal] == animals(a));
    %n = length(idcs);
    n = length(find(~isnan(data(1,idcs)) | ~isnan(data(2,idcs)))); % Cave, this assumes unpaired data - check on input side!
    colID = find(lut(:,1) == cellinfo(idcs(1)).animal);
    
    scat = scatter(nanmean(data(1,idcs)), nanmean(data(2,idcs)), .1+n*markerScaleFct, lut(colID,2:4),'filled'); % +.1 if n == 0
    scat.MarkerFaceAlpha = alpha;
    
    errorbar(nanmean(data(1,idcs)), nanmean(data(2,idcs)),...
        nanstd(data(2,idcs)),nanstd(data(2,idcs)),...
        nanstd(data(1,idcs)),nanstd(data(1,idcs)),...
        'CapSize',0,'LineWidth',LineWidth,'color','k');
    
    animal_means(a).animal = animals(a);
    animal_means(a).n = n;
    animal_means(a).data = nanmean(data(:,idcs),2);
    animal_means(a).SEM = nanstd(data(:,idcs),[],2);
end

end