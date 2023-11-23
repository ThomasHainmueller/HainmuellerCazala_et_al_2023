%% Utility functions for plotting data from the 'boutons dataset'
% To start, load all bouton data as CA1bt, CA3bt, DGbt
% Then load principal cell data as CA1pc, etc.

%% Make comparative cdfplots for boutons, pcs, fam and novel
% set e.g. property = 'AUCrate'; limits = [0,4]

vCA1pc = metastatistics(CA1pc,'cell',property,'rescategories',1:2,'categories',1:2);
vCA1bt = metastatistics(CA1bt,'cell',property,'rescategories',1:2,'categories',1:2);
vCA3pc = metastatistics(CA3pc,'cell',property,'rescategories',1:2,'categories',1:2);
vCA3bt = metastatistics(CA3bt,'cell',property,'rescategories',1:2,'categories',1:2);
vDGpc = metastatistics(DGpc,'cell',property,'rescategories',1:2,'categories',1:2);
vDGbt = metastatistics(DGbt,'cell',property,'rescategories',1:2,'categories',1:2);
close all;

figure;
subplot(3,1,1); hold on; xlim(limits);
h(1,1) = cdfplot(vCA1pc(:,1)); set(h(1,1),'LineStyle','-.','Color',[0,0,0]);
h(1,2) = cdfplot(vCA1pc(:,2)); set(h(1,2),'LineStyle','-.','Color',[.5,.5,.5]);
h(1,3) = cdfplot(vCA1bt(:,1)); set(h(1,3),'Color',[0,0,0]);
h(1,4) = cdfplot(vCA1bt(:,2)); set(h(1,4),'Color',[.5,.5,.5]);
title('CA1');

subplot(3,1,2); hold on; xlim(limits);
h(2,1) = cdfplot(vCA3pc(:,1)); set(h(2,1),'LineStyle','-.','Color',[0,0,1]);
h(2,2) = cdfplot(vCA3pc(:,2)); set(h(2,2),'LineStyle','-.','Color',[.5,.5,1]);
h(2,3) = cdfplot(vCA3bt(:,1)); set(h(2,3),'Color',[0,0,1]);
h(2,4) = cdfplot(vCA3bt(:,2)); set(h(2,4),'Color',[.5,.5,1]);
title('CA3'); 

subplot(3,1,3); hold on; xlim(limits);
h(3,1) = cdfplot(vDGpc(:,1)); set(h(3,1),'LineStyle','-.','Color',[1,0,0]);
h(3,2) = cdfplot(vDGpc(:,2)); set(h(3,2),'LineStyle','-.','Color',[1,.5,.5]);
h(3,3) = cdfplot(vDGbt(:,1)); set(h(3,3),'Color',[1,0,0]);
h(3,4) = cdfplot(vDGbt(:,2)); set(h(3,4),'Color',[1,.5,.5]);
title('DG'); 

clear vCA1pc vCA1bt vCA3pc vCA3pc vDGpc vDGbt;
%% Calculate spatial information divided by AUC 
for n=1:length(CA1)
    for c=1:length(CA1{n}.cells)
        CA1{n}.cells{c}.SIperAUC = ...
            CA1{n}.cells{c}.spatialinfo./CA1{n}.cells{c}.AUCrate;
    end
end

for n=1:length(CA3)
    for c=1:length(CA3{n}.cells)
        CA3{n}.cells{c}.SIperAUC = ...
            CA3{n}.cells{c}.spatialinfo./CA3{n}.cells{c}.AUCrate;
    end
end

for n=1:length(DG)
    for c=1:length(DG{n}.cells)
        DG{n}.cells{c}.SIperAUC = ...
            DG{n}.cells{c}.spatialinfo./DG{n}.cells{c}.AUCrate;
    end
end

clear n c
%% Plot a given property: 1. cdf, 2. boxplot
SI = 1;
PF = 1;
active = 0;

vCA1f = metastatistics(CA1,'cell',property,'rescategories',1,'categories',1,'plotting',0,...
    'active',active,'SI',SI,'PF',PF,'combination','or');
vCA1n = metastatistics(CA1,'cell',property,'rescategories',2,'categories',2,'plotting',0,...
    'active',active,'SI',SI,'PF',PF,'combination','or'); vCA1n = vCA1n(:,2);
vCA3f = metastatistics(CA3,'cell',property,'rescategories',1,'categories',1,'plotting',0,...
    'active',active,'SI',SI,'PF',PF,'combination','or');
vCA3n = metastatistics(CA3,'cell',property,'rescategories',2,'categories',2,'plotting',0,...
    'active',active,'SI',SI,'PF',PF,'combination','or'); vCA3n = vCA3n(:,2);
vDGf = metastatistics(DG,'cell',property,'rescategories',1,'categories',1,'plotting',0,...
    'active',.05,'combination','and');
vDGn = metastatistics(DG,'cell',property,'rescategories',2,'categories',2,'plotting',0,...
    'active',active,'SI',SI,'PF',PF,'combination','or'); vDGn = vDGn(:,2);

% Figure 1: cdf
figure; hold on; xlim(limits);
h(1,1) = cdfplot(vCA1f); set(h(1,1),'Color',[0,0,0]);
h(1,2) = cdfplot(vCA1n); set(h(1,2),'Color',[.5,.5,.5]);
h(1,3) = cdfplot(vCA3f); set(h(1,3),'Color',[0,0,1]);
h(1,4) = cdfplot(vCA3n); set(h(1,4),'Color',[.5,.5,1]);
h(1,5) = cdfplot(vDGf); set(h(1,5),'Color',[1,0,0]);
h(1,6) = cdfplot(vDGn); set(h(1,6),'Color',[1,.5,.5]);
title(property);

% Figure 2: boxplot
res = catuneven({vCA1f',vCA1n',vCA3f',vCA3n',vDGf',vDGn'},NaN);
figure; 
boxplot(res,'ColorGroup',[2,2,3,3,1,1],'Labels',{'CA1f','CA1n','CA3f','CA3n','DGf','DGn'});
title(property);

% Get p-values from rank-sum test
for r = 1:size(res,2)
    for c = 1:size(res,2)
        prank(r,c) = ranksum(res(:,r),res(:,c));
    end
end
clear vCA1f vCA1n vCA3f vCA3n vDGf vDGn h SI PF active

%% Experiment- and animal average
% Get animal#, x-offset and colorcode and  for each dataset.
SI = 1;
PF = 1;
active = 0;

for n=1:length(CA1)
    CA1av(n) = CA1{n}.animal;
    CA1xs(n) = mod(CA1av(n),10)/20-.25;
    CA1col(n,1:3) = [mod(CA1av(n),10)/10, mod(CA1av(n),7)/7, mod(CA1av(n),3)/3];
end
for n=1:length(CA3)
    CA3av(n) = CA3{n}.animal;
    CA3xs(n) = mod(CA3av(n),10)/20-.25;
    CA3col(n,1:3) = [mod(CA3av(n),10)/10, mod(CA3av(n),7)/7, mod(CA3av(n),3)/3];
end
for n=1:length(DG)
    DGav(n) = DG{n}.animal;
    DGxs(n) = mod(DGav(n),10)/20-.25;
    DGcol(n,1:3) = [mod(DGav(n),10)/10, mod(DGav(n),7)/7, mod(DGav(n),3)/3];    
end

% Get x-offset and colorcode for animal averages
CA1animals = unique(CA1av);
for n=1:length(CA1animals)
    CA1axs(n) = mod(CA1animals(n),10)/20-.25;
    CA1acol(n,1:3) = [mod(CA1animals(n),10)/10, mod(CA1animals(n),7)/7, mod(CA1animals(n),3)/3];
end

CA3animals = unique(CA3av);
for n=1:length(CA3animals)
    CA3axs(n) = mod(CA3animals(n),10)/20-.25;
    CA3acol(n,1:3) = [mod(CA3animals(n),10)/10, mod(CA3animals(n),7)/7, mod(CA3animals(n),3)/3];
end

DGanimals = unique(DGav);
for n=1:length(DGanimals)
    DGaxs(n) = mod(DGanimals(n),10)/20-.25;
    DGacol(n,1:3) = [mod(DGanimals(n),10)/10, mod(DGanimals(n),7)/7, mod(DGanimals(n),3)/3];
end

CA1anf = metastatistics(CA1,'animal',property,'categories',1,'rescategories',1,'plotting',0,...
    'active',active,'SI',SI,'PF',PF,'combination','or');
CA3anf = metastatistics(CA3,'animal',property,'categories',1,'rescategories',1,'plotting',0,...
    'active',active,'SI',SI,'PF',PF,'combination','or');
DGanf = metastatistics(DG,'animal',property,'categories',1,'rescategories',1,'plotting',0,...
    'active',active,'SI',SI,'PF',PF,'combination','or');
CA1ann = metastatistics(CA1,'animal',property,'categories',2,'rescategories',1:2,'plotting',0,...
    'active',active,'SI',SI,'PF',PF,'combination','or'); CA1ann = CA1ann(:,2);
CA3ann = metastatistics(CA3,'animal',property,'categories',2,'rescategories',1:2,'plotting',0,...
    'active',active,'SI',SI,'PF',PF,'combination','or'); CA3ann = CA3ann(:,2);
DGann = metastatistics(DG,'animal',property,'categories',2,'rescategories',1:2,'plotting',0,...
    'active',active,'SI',SI,'PF',PF,'combination','or'); DGann = DGann(:,2);
anmean = catuneven({CA1anf',CA1ann',CA3anf',CA3ann',DGanf',DGann'},NaN);

CA1dsf = metastatistics(CA1,'hybrid',property,'categories',1,'rescategories',1,'plotting',0,...
    'active',active,'SI',SI,'PF',PF,'combination','or');
CA3dsf = metastatistics(CA3,'hybrid',property,'categories',1,'rescategories',1,'plotting',0,...
    'active',active,'SI',SI,'PF',PF,'combination','or');
DGdsf = metastatistics(DG,'hybrid',property,'categories',1,'rescategories',1,'plotting',0,...
    'active',active,'SI',SI,'PF',PF,'combination','or');
CA1dsn = metastatistics(CA1,'hybrid',property,'categories',2,'rescategories',1:2,'plotting',0,...
    'active',active,'SI',SI,'PF',PF,'combination','or'); CA1dsn = CA1dsn(:,2);
CA3dsn = metastatistics(CA3,'hybrid',property,'categories',2,'rescategories',1:2,'plotting',0,...
    'active',active,'SI',SI,'PF',PF,'combination','or'); CA3dsn = CA3dsn(:,2);
DGdsn = metastatistics(DG,'hybrid',property,'categories',2,'rescategories',1:2,'plotting',0,...
    'active',active,'SI',SI,'PF',PF,'combination','or'); DGdsn = DGdsn(:,2);
% CA1ds = metastatistics(CA1,'hybrid',property,'categories',1:2,'rescategories',1:2,'plotting',0);
% CA3ds = metastatistics(CA3,'hybrid',property,'categories',1:2,'rescategories',1:2,'plotting',0);
% DGds = metastatistics(DG,'hybrid',property,'categories',1:2,'rescategories',1:2,'plotting',0);
dsmean = catuneven({CA1dsf',CA1dsn',CA3dsf',CA3dsn',DGdsf',DGdsn'},NaN);

% Toggle whether familiar and novel are displayed by creating a variable named
% 'both' in the workspace.
if ~exist('both')
    figure; hold on;
    % Plot dataset means
    scatter([ones(1,length(CA1dsf))*1.+CA1xs],CA1dsf,[],CA1col,'x');
    scatter([ones(1,length(CA3dsf))*2.+CA3xs],CA3dsf,[],CA3col,'x');
    scatter([ones(1,length(DGdsf))*3.+DGxs],DGdsf,[],DGcol,'x');
    
    % Plot animal means
    scatter([ones(1,length(CA1animals))*1.+CA1axs],CA1anf,[],CA1acol,'filled');
    scatter([ones(1,length(CA3animals))*2.+CA3axs],CA3anf,[],CA3acol,'filled');
    scatter([ones(1,length(DGanimals))*3.+DGaxs],DGanf,[],DGacol,'filled');
    
    title(property);
    xlabel(['CA1 -- CA3 -- DG']);
else
    figure; hold on;
    % Plot dataset means
    for n=1:length(CA1dsf)
        plot([1+CA1xs(n) 2+CA1xs(n)],[CA1dsf(n) CA1dsn(n)],'Color',CA1col(n,:),'Marker','.');
    end
    for n=1:length(CA3dsf)
        plot([3+CA3xs(n) 4+CA3xs(n)],[CA3dsf(n) CA3dsn(n)],'Color',CA3col(n,:),'Marker','.');
    end
    for n=1:length(DGdsf)
        plot([5+DGxs(n) 6+DGxs(n)],[DGdsf(n) DGdsn(n)],'Color',DGcol(n,:),'Marker','.');
    end
    
    % Plot animal means
    for n=1:length(CA1anf)
        plot([1+CA1axs(n) 2+CA1axs(n)],[CA1anf(n) CA1ann(n)],'Color',CA1acol(n,:),'Marker','o');
    end
    for n=1:length(CA3anf)
        plot([3+CA3axs(n) 4+CA3axs(n)],[CA3anf(n) CA3ann(n)],'Color',CA3acol(n,:),'Marker','o');
    end
    for n=1:length(DGanf)
        plot([5+DGaxs(n) 6+DGaxs(n)],[DGanf(n) DGann(n)],'Color',DGacol(n,:),'Marker','o');
    end
    
    title(property);
    xlabel(['CA1f - CA1n - CA3f - CA3n - DGf - DGn']);
end

clear CA1an CA1ds CA1xs CA1col CA1axs CA1acol CA3an CA3ds CA3xs CA3col CA3axs...
    CA3acol DGan DGds DGxs DGcol DGaxs DGacol CA1animals CA3animals DGanimals...
    n CA1av CA3av DGav CA1anf CA1ann CA3anf CA3ann DGanf DGann...
    CA1dsf CA1dsn CA3dsf CA3dsn DGdsf DGdsn active SI PF
    
%% Split datasets to blocks and calculate correlations.

for n=1:length(CA1)
    CA1s{n} = block_categories(CA1{n},1:2,5);
    CA1s{n} = spatial_corr(CA1s{n},1,0,1,1);
end

for n=1:length(CA3)
    CA3s{n} = block_categories(CA3{n},1:2,5);
    CA3s{n} = spatial_corr(CA3s{n},1,0,1,1);
end

for n=1:length(DG)
    DGs{n} = block_categories(DG{n},1:2,5);
    DGs{n} = spatial_corr(DGs{n},1,0,1,1);
end

%% Make cross-correlation boxplots
CA1corr = metacorr(CA1s,'SI',.05,'active',.05,'categories',1:6,'combination','or');
CA3corr = metacorr(CA3s,'SI',.05,'active',.05,'categories',1:6,'combination','or');
DGcorr = metacorr(DGs,'SI',.05,'active',.05,'categories',1:6,'combination','or');

allcorr = catuneven({CA1corr(:,1,3)',CA1corr(:,2,4)',CA1corr(:,1,2)',...
    CA3corr(:,1,3)',CA3corr(:,2,4)',CA3corr(:,1,2)',DGcorr(:,1,3)',...
    DGcorr(:,2,4)',DGcorr(:,1,2)'},NaN);

figure; boxplot(allcorr,'ColorGroup',[2,2,2,3,3,3,1,1,1],'Labels',{...
    'CA1-FF','CA1-NN', 'CA1-FN', 'CA3-FF','CA3-NN','CA3-FN',...
    'DG-FF','DG-NN','DG-FN'});
clear CA1corr CA3corr DGcorr
%% Calculate mean speed per track
allspeeds = [];

for n=1:length(CA1)
    for c=1:length(CA1{n}.metadata.categories)
        for r = 1:length(CA1{n}.metadata.categories{c}.y)
            % Multiply by acquisition rate and 200 (2V/400cm) to get cm/s
            thisspeed = abs(diff(CA1{n}.metadata.categories{c}.y{r}))...
                .*CA1{n}.metadata.categories{c}.acquisition_rate(r).*200;
            % Correct for excessive 'jump' values
            thisspeed(thisspeed>prctile(thisspeed,98))=0;
            speeds(c,r) = nanmean(thisspeed);
        end
    end
    CA1{n}.metadata.speed = nanmean(speeds,2);
    allspeeds(end+1,:) = CA1{n}.metadata.speed;
end

for n=1:length(CA3)
    for c=1:length(CA3{n}.metadata.categories)
        for r = 1:length(CA3{n}.metadata.categories{c}.y)
            % Multiply by acquisition rate and 200 (2V/400cm) to get cm/s
            thisspeed = abs(diff(CA3{n}.metadata.categories{c}.y{r}))...
                .*CA3{n}.metadata.categories{c}.acquisition_rate(r).*200;
            % Correct for excessive 'jump' values
            thisspeed(thisspeed>prctile(thisspeed,98))=0;
            speeds(c,r) = nanmean(thisspeed);
        end
    end
    CA3{n}.metadata.speed = nanmean(speeds,2);
    allspeeds(end+1,:) = CA3{n}.metadata.speed;
end

for n=1:length(DG)
    for c=1:length(DG{n}.metadata.categories)
        for r = 1:length(DG{n}.metadata.categories{c}.y)
            % Multiply by acquisition rate and 200 (2V/400cm) to get cm/s
            thisspeed = abs(diff(DG{n}.metadata.categories{c}.y{r}))...
                .*DG{n}.metadata.categories{c}.acquisition_rate(r).*200;
            % Correct for excessive 'jump' values
            thisspeed(thisspeed>prctile(thisspeed,98))=0;
            speeds(c,r) = nanmean(thisspeed);
        end
    end
    DG{n}.metadata.speed = nanmean(speeds,2);
    allspeeds(end+1,:) = DG{n}.metadata.speed;
end

figure; plot(allspeeds','Marker','o');
hold on; errorbar(nanmean(allspeeds,1),nanstd(allspeeds,[],1)./sqrt(size(allspeeds,1)));

clear n r c thisspeed speeds 
%% Get numbers of active- and spatially modulated cells
for n=1:length(CA1)
    ncells = length(CA1{n}.cells);
    nactive = length(findcells(CA1{n},1,.05,1,1));
    nSI = length(findcells(CA1{n},1,0.05,.05,1));
    nPF = length(findcells(CA1{n},1,0,1,.05));
    FactCA1(n) = nactive/ncells;
    FsiCA1(n) = nSI/ncells;
    FpfCA1(n) = nPF/ncells;
end

for n=1:length(CA3)
    ncells = length(CA3{n}.cells);
    nactive = length(findcells(CA3{n},1,.05,1,1));
    nSI = length(findcells(CA3{n},1,.05,.05,1));
    nPF = length(findcells(CA3{n},1,0,1,.05));
    FactCA3(n) = nactive/ncells;
    FsiCA3(n) = nSI/ncells;
    FpfCA3(n) = nPF/ncells;
end

for n=1:length(DG)
    ncells = length(DG{n}.cells);
    nactive = length(findcells(DG{n},1,.05,1,1));
    nSI = length(findcells(DG{n},1,.05,.05,1));
    nPF = length(findcells(DG{n},1,0,1,.05));
    FactDG(n) = nactive/ncells;
    FsiDG(n) = nSI/ncells;
    FpfDG(n) = nPF/ncells;
end

res = catuneven({FactCA1,FactCA3,FactDG,FsiCA1,FsiCA3,FsiDG,...
    FpfCA1,FpfCA3,FpfDG},NaN);
figure; hold on;
bar(nanmean(res,1));
% CAVE: Lenght of datasets hard-coded for convenience (SEM)
errorbar(nanmean(res,1),nanstd(res,[],1)./sqrt([9 11 9 9 11 9 9 11 9]),'x');
plot(res','o','Color',[.5 .5 .5]);
title('Fraction active - Fraction active&SI - Fraction PF');

clear ncells nactive nSI nPF n
%% Plot 5 random place cells and significant SI cells
for n=1:5
    CA1plc = findcells(CA1{n},1,0,1,.05);
    CA1sic = findcells(CA1{n},1,.05,.05,1);
    CA3plc = findcells(CA3{n},1,0,1,.05);
    CA3sic = findcells(CA3{n},1,.05,.05,1);
    DGplc = findcells(DG{n},1,0,1,.05);
    DGsic = findcells(DG{n},1,.05,.05,1);
    
    plot_placeactivity(CA1{n},CA1plc(randi(length(CA1plc))),1);
    title('CA1 place cell');
    plot_placeactivity(CA3{n},CA3plc(randi(length(CA3plc))),1);
    title('CA3 place cell');
    plot_placeactivity(DG{n},DGplc(randi(length(DGplc))),1);
    title('DG place cell');
    
    plot_placeactivity(CA1{n},CA1sic(randi(length(CA1sic))),1);
    title('CA1 SI cell');
    plot_placeactivity(CA3{n},CA3sic(randi(length(CA3sic))),1);
    title('CA3 SI cell');
    plot_placeactivity(DG{n},DGsic(randi(length(DGsic))),1);
    title('DG SI cell');
    end

clear CA1plc CA1sic CA3plc CA3sic DGplc DGsic n