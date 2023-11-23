function data = spatial_param_rnd(data,varargin)
% Calculate [Moving/Immobile ratio, speed-modulation slope], SI, SI/AUC, 
% [trial-trial variability], 1st to 2nd half stability, tuning-vector length 
% for real- and circularly shuffled data.

p = inputParser;
addParameter(p,'nshuffles',1000,@isnumeric);
addParameter(p,'bins',.1:.025:1.9,@isnumeric);
addParameter(p,'category',1,@isnumeric);

parse(p,varargin{:})

nshuffles = p.Results.nshuffles;
bins = p.Results.bins;
ca = p.Results.category;

rng(1);
framerate = data.metadata.categories{1}.acquisition_rate(1);

%% Generate LUT for extracting per-bin information
mvg = logical(cat(2,data.metadata.categories{ca}.moving{:}));
y = cat(2,data.metadata.categories{ca}.y{:});
y = y(mvg);

% Create circular coordinates for tuning vector length (TVL)
xrad = ((1:length(bins)-1)*2*pi/(length(bins)-1))';

% Create logical matrix to discretize position
for b = length(bins):-1:2
    yb(b-1,:) = y>bins(b-1) & y<=bins(b);
end

ryb = false(size(yb,1),size(yb,2),nshuffles);
ry = zeros(length(y),nshuffles);

parfor sh = 1:nshuffles
    ry(:,sh) = circshift(y,round(rand*length(y)));
end

%% Calculate random spatial maps and extract information
for n = 1:length(data.cells)
    tic
    fprintf('Processing cell# %d\n',n);
    tr = cat(2,data.cells{n}.categories{ca}.dFoT{:});
    tr = tr(mvg);
    
    % Remove NaN from trace - Added 221129
    tr(isnan(tr)) = 0;
    
    mid = round(size(tr,2)/2);
    
    dFoY = discretize(tr,y,bins);
    dFoY1 = discretize(tr(1:mid),y(1:mid),bins);
    dFoY2 = discretize(tr(mid:end),y(mid:end),bins);
    
    % Interpolate missing bins, if present - introduced 22/11/29
    dFoY = interp1(bins(~isnan(dFoY)),dFoY(~isnan(dFoY)),bins(1:end-1),'linear',nanmean(dFoY))';
    dFoY1 = interp1(bins(~isnan(dFoY1)),dFoY1(~isnan(dFoY1)),bins(1:end-1),'linear',nanmean(dFoY1))';
    dFoY2 = interp1(bins(~isnan(dFoY2)),dFoY2(~isnan(dFoY2)),bins(1:end-1),'linear',nanmean(dFoY2))';
    
    % Spatial info
    SI = spatial_info(tr,y,bins,framerate);
    SIpAUC = SI/data.cells{n}.AUCrate(ca);
    
    % Spatial coherence 
    Coh = calc_coherence(dFoY);
    
    % Stability
    Stab = corrcoef(dFoY1,dFoY2);
    Stab = Stab(1,2);
    
    % Circular tuning vector
    TVL = circ_r(xrad,dFoY);
    TVA = circ_mean(xrad,dFoY);

%     figure;
%     polarplot(xrad,dFoY);
%     hold on;
%     polarscatter(TVA,TVL);

    parfor sh = 1:nshuffles
        dFoYr = discretize(tr,ry(:,sh),bins);
        dFoYr1 = discretize(tr(1:mid),ry(1:mid,sh),bins);
        dFoYr2 = discretize(tr(mid:end),ry(mid:end,sh),bins);
        
        % Interpolate missing bins, if present - introduced 22/11/29
        dFoYr = interp1(bins(~isnan(dFoYr)),dFoYr(~isnan(dFoYr)),...
            bins(1:end-1),'linear',nanmean(dFoYr))';
        dFoYr1 = interp1(bins(~isnan(dFoYr1)),dFoYr1(~isnan(dFoYr1)),...
            bins(1:end-1),'linear',nanmean(dFoYr1))';
        dFoYr2 = interp1(bins(~isnan(dFoYr2)),dFoYr2(~isnan(dFoYr2)),...
            bins(1:end-1),'linear',nanmean(dFoYr2))';
        
        SIr(sh) = spatial_info(tr,ry(:,sh),bins,framerate);
        SIpAUC(sh) = SIr(sh)/data.cells{n}.AUCrate(ca);
        thisStab = corrcoef(dFoYr1,dFoYr2);
        Stabr(sh) = thisStab(1,2);
        Cohr(sh) = calc_coherence(dFoYr);
        TVLr(sh) = circ_r(xrad,dFoYr);
    end
    
    % Calculate p-values
    pSI = length(find(SIr >= SI))/nshuffles;
    pStab = length(find(Stabr >= Stab))/nshuffles;
    pCoh = length(find(Cohr >= Coh))/nshuffles;
    pTVL = length(find(TVLr >= TVL))/nshuffles;
    
    % Store results, random value(s) and p-values with data.
    data.cells{n}.spatial_P(ca) = pSI;
    data.cells{n}.spatialinfo(ca) = SI;
    data.cells{n}.spatialinfo_rand{ca} = SIr;
    
    data.cells{n}.SIperAUC = SI/data.cells{n}.AUCrate(ca); % ADDED 210301
    
    data.cells{n}.spatial_coherence(ca) = Coh;
    data.cells{n}.spatial_coherence_P(ca) = pCoh;
    data.cells{n}.spatial_coherence_rand{ca} = Cohr;
    
    data.cells{n}.tvl(ca) = TVL;
    data.cells{n}.tva(ca) = TVA;
    data.cells{n}.tvl_P(ca) = pTVL;
    data.cells{n}.TVLr{ca} = TVLr;
    
    data.cells{n}.sessionstab(ca) = Stab;
    data.cells{n}.sessionstab_P(ca) = pStab;
    data.cells{n}.sessionstab_rand{ca} = Stabr;

    % Adapt these for plotting vs. random distribution - cave, sets only
    % one value derived from the current category! - ADDED 221129
    data.cells{n}.spatialinfo_r(1) = SIr(1);
    data.cells{n}.SIperAUC_r(1) = SIr(1) / data.cells{n}.AUCrate(ca);
    data.cells{n}.spatial_coherence_r(1) = Cohr(1);
    data.cells{n}.sessionstab_r(1) = Stabr(1);
    data.cells{n}.tvl_r(1) = TVLr(1);

    
%     figure; imagesc(dFoYr);
%     hold on; plot(dFoY*1000,'k');
%     plot(nanmean(dFoYr,1)*10000);
% 
%     figure; hold on;
%     plot(nanmean(dFoYr,1),'color',[.5 .5 .5]);
%     plot(dFoY,'k'); 
%     plot(dFoY1,'r');
%     plot(dFoY2,'b');
    toc
end
end

function Coh = calc_coherence(dFoY)

for f = length(dFoY)-1:-1:2
    % Get values surrounding each position on the trace
    dFoYsurr(f-1) = nanmean([dFoY(f+1),dFoY(f-1)]);
end

c = corrcoef(dFoY(2:end-1),dFoYsurr);
Coh = c(1,2);
end