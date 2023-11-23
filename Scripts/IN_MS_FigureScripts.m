%% REPLICATING FIGURES IN HAINMUELLER, CAZALA ET AL., 2023
% Requires the following data files to run (~19GB total): 
% 230318_CA1_DREADD, 230825_ClozapineControl_DG,
% 230901_DREADD_IN_DG, 230905_allIN, AxonData_TH_221206,
% DREADDdata_TH_231121

% FIGURE POSITION: [left bottom width height]


%% FIGURE 1
%% E-K) Representative examples with traces, map and speed modulation
names = {'CA1PV','CA1SOM','CA3PV','CA3SOM','DGPV','DGSOM'};
signals = {'dFoT','dFodY','dFoY'};%{'zscored','dZodY','dZoY'};
sig_fact = [2 4]; % Scaling factor for traces, 0.5 means 0.5 DF/F per plotted unit
ex = [1 6 2; 6 4 3];
VarSaveStrings = {'Fig1f','Fig1j'};
xlims = [250 450; 150 350]; % Set time axis to respective windows
colmap = jet;
colmap_lim = [0 1.5; 0 .5];
f1 = figure;
set(f1,'defaultAxesFontSize',8,'Units','centimeters','Position',[2 2 13.2 8.2]);

for e = 1:size(ex,1)
    eval(sprintf('mds = %s;',names{ex(e,1)}))
    
    % SUBPLOT 1: TRACES OVER TIME
    subplot(size(ex,1),3,3*(e-1)+1,'Position',[.02 1.05-e/size(ex,1) .3 .9/size(ex,1)]); % [left bottom width height]
    hold on;
    fr = mds{ex(e,2)}.metadata.categories{1}.acquisition_rate(1);
    
    % Put together all traces
    [tr,v,y,dZoY] = deal([]);
    famI = 0;
    
    for r = 1:length(mds{ex(e,2)}.metadata.categories{1}.y)
        famL = length(mds{ex(e,2)}.cells{ex(e,3)}.categories{1}.(signals{1}){r});
        dtr = cat(2,mds{ex(e,2)}.cells{ex(e,3)}.categories{1}.(signals{1}){r},...
            mds{ex(e,2)}.cells{ex(e,3)}.categories{2}.(signals{1}){r});
        tr = cat(2,tr,dtr);
        v = cat(2,v,mds{ex(e,2)}.metadata.categories{1}.v{r},...
            mds{ex(e,2)}.metadata.categories{2}.v{r});
        y = cat(2,y,mds{ex(e,2)}.metadata.categories{1}.y{r},...
            mds{ex(e,2)}.metadata.categories{2}.y{r});
        
        % Plot rectangle for new track
        rectangle('Position',[famI/fr,-10,famL/fr,15],'FaceColor',[.7 .7 .7],...
            'LineWidth',1,'EdgeColor',[.7 .7 .7])
        
        famI = famI + length(dtr);
    end
    clear famI
    
    t = (0:length(v)-1)/fr;
    v(v>70) = 70; % Speed-scale: max = 30 cm/s
    tr = movmean(tr,10);
    
    % Calcium trace (black)
    plot(t,tr*sig_fact(e),'k');
    % Speed (blue)
    plot(t,v/20-2,'b');
    % Track location (y; red)
    plot(t,y-4,'r');
    
    %xlim([0 180]);
    xlim(xlims(e,:));
    xlabel('Time (s)');
    ylim([-4.5 4.5])
    
    % SUBPLOT 2: SPEED TUNING
    subplot(size(ex,1),3,3*(e-1)+2,'Position',[.4 1.08-e/size(ex,1) .25 .82/size(ex,1)]);
    hold on
    xv = 0:15;
    pro = nanmean(cat(2,mds{ex(e,2)}.cells{ex(e,3)}.categories{1}.(signals{2}){:}),2);
    prosem = nanstd(cat(2,mds{ex(e,2)}.cells{ex(e,3)}.categories{1}.(signals{2}){:}),[],2)./...
        sqrt(length(mds{ex(e,2)}.cells{ex(e,3)}.categories{1}.(signals{2})));
    sR = mds{ex(e,2)}.cells{ex(e,3)}.speed_R(1);
    sM = mds{ex(e,2)}.cells{ex(e,3)}.speedmod(1);
    sp = mds{ex(e,2)}.cells{ex(e,3)}.speed_P(1);
    errorbar(xv,pro(1:16),prosem(1:16),'k');
    P = polyfit(xv(~isnan(pro(1:16))),pro(~isnan(pro(1:16)))',1);
    %title(sprintf('M=%.2d, p=%.2d',sR,sp),'FontSize',6);
    plot(xv,xv*P(1)+P(2),'.','color','b');
    ylabel('Activity (DF/F)','FontSize',8);
    xlabel('Speed (cm/s)');
    
    % Save source data files
    T = table(xv',pro(1:16),prosem(1:16),'VariableNames',{'Speed (cm/s)','Mean Activity','SEM Activity'});
    writetable(T,VarSaveStrings{e},'Delimiter','\t');
    
    % SUBPLOT 3: SPATIAL MAPS
    subplot(size(ex,1),3,3*(e-1)+3,'Position',[.7 1.08-e/size(ex,1) .28 .87/size(ex,1)]);
    dZoY = cat(2,dZoY,mds{ex(e,2)}.cells{ex(e,3)}.categories{1}.(signals{3}){:});
    dZoY = movmean(dZoY,3,1);
    xv = 0:4;
    yv = 1:size(dZoY,2);
    
    imagesc(xv,yv,dZoY');
    caxis(colmap_lim(e,:));
    colormap(colmap)
    set(gca,'box','off');
    xlabel('Track distance (m)');
    %ylabel('Run #');
    
end

pos = get(f1,'Position');
set(f1,'PaperPositionMode','Auto','PaperUnits','Inches',...
    'PaperSize',[pos(3), pos(4)],'Renderer','Painters')
print(f1,'ExampleCells','-dpdf','-r0')

clear names v t tr pro P xv ax mds e ex prosem famI famL fr dtr r y yv...
    signals colmap colmap_lim xlims sig_fact pos sM sp sR dZoY


%% FIGURE 2
%% A) Speedmaps colored by significance
% See IN_scripts - Extract speed profiles
names = {'CA1PV','CA1SOM','CA3PV','CA3SOM','DGPV','DGSOM'};
labels = {'CA1-PV','CA1-SOM','CA2/3-PV','CA2/3-SOM','DG-PV','DG-SOM'};

colors = {[0 1 1],[1 1 0],[.5 .5 .5]};
%plo = [1 4 2 5 3 6];
plo = 1:6;
xv = (0:15)';

%Make color bar figure
colbar(1,:,:) = ones(256,3)-repmat(0:1/255:1,3,1)'.*colors{1};
colbar(2,:,:) = ones(256,3)-repmat(0:1/255:1,3,1)'.*colors{2};
colbar(3,:,:) = ones(256,3)-repmat(0:1/255:1,3,1)'.*colors{3};
figure; image(colbar);
box off

ff = figure; sgtitle('DF/F');
set(ff,'Units','Inches','defaultAxesFontSize',8,'Position',[2 2 5 5.6]);
fz = figure; sgtitle('z-scored');
set(fz,'Units','Inches','defaultAxesFontSize',8,'Position',[6 2 4.5 5.6]);
fzc = figure; sgtitle('z-scored & scaled');
set(fzc,'Units','Inches','defaultAxesFontSize',8,'Position',[10 2 3.5 5.6]);
ffc = figure; sgtitle('DF/F & scaled');
set(ffc,'Units','Inches','defaultAxesFontSize',8,'Position',[10 2 3.5 5.6]);

for nmds = 1:length(names)
    eval(sprintf('mds = %s;',names{nmds}))
    gn = 0;
    for nds = 1:length(mds)
        for n = 1:length(mds{nds}.cells)
            gn = gn+1;
            for c = 1:length(mds{nds}.metadata.categories)
                ftr = nanmean(cat(2,...
                    mds{nds}.cells{n}.categories{c}.dFodY{:}),2);
                ztr = nanmean(cat(2,...
                    mds{nds}.cells{n}.categories{c}.dZodY{:}),2);
                ftr = ftr(1:16); ztr = ztr(1:16);
                
                ftr = movmean(ftr,3,'omitnan');
                ztr = movmean(ztr,3,'omitnan');
                fnan = isnan(ftr); znan = isnan(ztr);
                
                DF(:,gn,c) = ftr;
                DZ(:,gn,c) = ztr;
                
                Pfit = polyfit(xv(~fnan),ftr(~fnan),1);
                mds{nds}.cells{n}.speedmod(c) = Pfit(1);
                mds{nds}.cells{n}.speedYincept(c) = Pfit(2);
                
                % Correlation
                R = mds{nds}.cells{n}.speed_R(c);
                p = mds{nds}.cells{n}.speed_P(c);
                
                % Make colored plot
                ztr = ztr-min(ztr);
                ztr = ztr/max(ztr);
                ftr = ftr-min(ftr);
                ftr = ftr/max(ftr);
                
                
                if p < .05 && R > 0
                    DZc(gn,:,:,c) = ones(size(ztr,1),3)-ztr*colors{1};
                    DFc(gn,:,:,c) = ones(size(ftr,1),3)-ftr*colors{1};
                elseif p < .05 && R < 0
                    DZc(gn,:,:,c) = ones(size(ztr,1),3)-ztr*colors{2};
                    DFc(gn,:,:,c) = ones(size(ftr,1),3)-ftr*colors{2};
                else
                    DZc(gn,:,:,c) = ones(size(ztr,1),3)-ztr*colors{3};
                    DFc(gn,:,:,c) = ones(size(ftr,1),3)-ftr*colors{3};
                end
                
                %order{nmds}(gn,c) = R(1,2);
                order{nmds}(gn,c) = Pfit(1);
            end
        end
    end
    
    for c = 1:length(mds{nds}.metadata.categories)
        [~,order{nmds}(:,c)] = sort(order{nmds}(:,c));
    end
    
    DFplot{nmds} = DF;
    DZplot{nmds} = DZ;
    DZcplot{nmds} = DZc;
    DFcplot{nmds} = DFc;
    
    yv = (1:size(DFplot{nmds},2))';
    
    figure(fz);
    subplot(3,2,plo(nmds));
    imagesc(xv,yv,DZplot{nmds}(:,order{nmds}(:,1),1)'); 
    caxis([-.5 1.5]); colormap(hot); xlim([0 15]);
    xlabel('Speed (cm/s)'); ylabel('Cell #');
    title(names{nmds});
    colorbar;
    
    figure(ff);
    subplot(3,2,plo(nmds));
    imagesc(xv,yv,DFplot{nmds}(:,order{nmds}(:,1),1)'); 
    caxis([0 .7]); xlim([0 15]);
    xlabel('Speed (cm/s)'); ylabel('Cell #');
    title(names{nmds});
    colorbar
    
    figure(fzc);
    subplot(3,2,plo(nmds));
    image(xv,yv,DZcplot{nmds}(order{nmds}(:,1),:,:,1)); 
    xlabel('Speed (cm/s)'); ylabel('Cell #');
    title(names{nmds});
    yticklabels('');
    unitstr = sprintf('N = %i',gn);
    text(0.3,gn-gn/11,unitstr)
    ylabel('Cell #');
    title(labels{nmds});
    
    figure(ffc);
    subplot(3,2,plo(nmds));
    image(xv,yv,DFcplot{nmds}(order{nmds}(:,1),:,:,1)); 
    xlabel('Speed (cm/s)'); ylabel('Cell #');
    title(names{nmds});
    yticklabels('');
    unitstr = sprintf('N = %i',gn);
    text(0.3,gn-gn/11,unitstr)
    ylabel('Cell #');
    title(labels{nmds});
    
    % WRITE BACK TO ORIGINAL DATA (.speedmod and .speedYintercep are added
    % to cell).
    %eval(sprintf('%s = mds;',names{nmds}))
    
    clear DF DZ DZc DFc
end
clear mds nds nmds n c r names DF DZ DFc gn fi fnan ftr P xv yv znan ztr ...
    colors DFcplot DFplot DZplot f1 f2 f3 p plo R

%% B) Speed tuning pie charts and Chi2test
names = {'CA1PV','CA1SOM','CA3PV','CA3SOM','DGPV','DGSOM'};
labels = {'CA1-PV','CA1-SOM','CA2/3-PV','CA2/3-SOM','DG-PV','DG-SOM'};
f1 = figure;%('Position',[200,100,120,800]);
set(f1,'Units','Inches','Position',[2 1 1.2 6],'defaultAxesFontSize',12);

for nmds = 1:length(names)
    eval(sprintf('mds = %s;',names{nmds}))
    speed_P = metastatistics(mds,'cell','speed_P','categories',1,...
        'rescategories',1,'plotting',false)';
    speed_R = metastatistics(mds,'cell','speed_R','categories',1,...
        'rescategories',1,'plotting',false)';
    tuned_Pos = length(find(speed_P<.05 & speed_R>0));
    tuned_Neg = length(find(speed_P<.05 & speed_R<0));
    untuned = length(find(speed_P>.05));
    
    tbl(1:3,nmds) = cat(1,tuned_Pos, tuned_Neg, untuned);
    plbs = tbl(1:3,nmds)/sum(tbl(1:3,nmds))*100;
    plbs = {sprintf('%.1f%%',plbs(1)), sprintf('%.1f%%',plbs(2)),...
        sprintf('%.1f%%',plbs(3))};
    
    figure(f1);
    subplot(6,1,nmds)
    pie(tbl(:,nmds),plbs);
    ax = gca;
    ax.Colormap = [1 0 0; 0 0 1; .5 .5 .5];
    title(labels{nmds},'FontSize',8);
    clear plbs
end

set(findobj(f1,'type','text'),'fontsize',7)

% Save source data files
PCtbl = tbl(1:3,:)./sum(tbl(1:3,:))*100;
T = table(PCtbl(1,:)',PCtbl(2,:)',PCtbl(3,:)',...
    'VariableNames',{'Pos%','Neg%','None%'},'RowNames',names);
writetable(T,'Fig2b','Delimiter','\t','WriteRowNames',true)  

% Save as pdf
pos = get(f1,'Position');
set(f1,'PaperPositionMode','Auto','PaperUnits','Inches',...
    'PaperSize',[pos(3), pos(4)],'Renderer','Painters')
print(f1,'SpeedPiecharts.pdf','-dpdf','-r0')

[Chi2,p] = chi2tbltest(tbl);
fprintf('Chi2 = %.2f; p = %d\n',Chi2,p);
clear names speed_P speed_R tuned_Pos tuned_Neg untuned tbl nmds mds ax Chi2 p f1 labels

%% C) Speed modulation violin plots
names = {'CA1PV','CA1SOM','CA3PV','CA3SOM','DGPV','DGSOM'};
param = 'speedmod';%'MovImmRatio','speedmod','trialvar','sessionstab','spatial_coherence','ZdiffScore'
SI = 1;
yl = 'Speed modulation slope';%param; 'Activity Ratio moving/immobile';
OutlierThr = 3; % In interquartile ranges, use tukey procedure
colors = {[0 .5 0],[.5 1 .5],[0 0 1],[.5 .5 1],[1 0 0],[1 .5 .5]};

for nmds = 1:length(names)
    eval(sprintf('mds = %s;',names{nmds}))
    val = metastatistics(mds,'cell',param,'categories',1,...
        'rescategories',1,'plotting',false,'SI',SI,'plotting',false)';
    % Use tukey outlier removal
    val = tukeyOutlierRemoval(val',OutlierThr);
    values{nmds} = val;
end
close all

f3 = figure; hold on;
set(f3,'Units','centimeters','position',[5 5 6.7 5.8],'defaultAxesFontSize',8);
plot([0 4],[0 0],':','color',[.7 .7 .7],'LineWidth',1.5);
beeswarmplot(catuneven(values(1:2:end),NaN),'orientation','left','color',colors(1:2:end),'markersize',6);
distributionPlot(values(1:2:end),'histOri','left','color',colors(1:2:end),'widthDiv',[2 1],'showMM',6,'FaceAlpha',.4);
beeswarmplot(catuneven(values(2:2:end),NaN),'orientation','right','color',colors(2:2:end),'markersize',6);
distributionPlot(values(2:2:end),'histOri','right','color',colors(2:2:end),'widthDiv',[2 2],'showMM',6,'FaceAlpha',.4);
xticklabels({'CA1','CA3','DG'});
ylabel(yl);
ylim([-.02 .04]);
yticks([-.02 0 .02 .04])

% Save source data files
VAL = catuneven(values,NaN);
T = table(VAL(:,1),VAL(:,2),VAL(:,3),VAL(:,4),VAL(:,5),VAL(:,6),...
    'VariableNames',names);
writetable(T,'Fig2c','Delimiter','\t','WriteRowNames',false)  

pos = get(f3,'Position');
set(f3,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(f3,sprintf('%s-SI%.2f.pdf',param,SI),'-dpdf','-r0')

kruskalwallisNdunns(catuneven(values,NaN),names);

clear param nmds mds yl names val Foutliers colors OutlierThr VAL

%% D) Activity Ratio moving/immobile
names = {'CA1PV','CA1SOM','CA3PV','CA3SOM','DGPV','DGSOM'};
param = 'MovImmRatio';%'MovImmRatio','speedmod','trialvar','sessionstab','spatial_coherence','ZdiffScore'
SI = 1;
yl = 'Log2 Activity Ratio mov/imm';
OutlierThr = 3; % In interquartile ranges, use tukey procedure
colors = {[0 .5 0],[.5 1 .5],[0 0 1],[.5 .5 1],[1 0 0],[1 .5 .5]};
%ylims = [.1 10];
ylims = [-3 3];

for nmds = 1:length(names)
    eval(sprintf('mds = %s;',names{nmds}))
    val = metastatistics(mds,'cell',param,'categories',1,...
        'rescategories',1,'plotting',false,'SI',SI,'plotting',false)';
    % Use tukey outlier removal
    val = tukeyOutlierRemoval(val',OutlierThr);
    
    KWvalues{nmds} = val;
    values{nmds} = real(log2(val));
    
    MvIm(:,1) = metastatistics(mds,'cell','Amov','plotting',false,'categories',1,'rescategories',1);
    MvIm(:,2) = metastatistics(mds,'cell','Aimm','plotting',false,'categories',1,'rescategories',1);
    MvIm = tukeyOutlierRemoval(MvIm,OutlierThr);
    p(nmds) = signrank(MvIm(:,1),MvIm(:,2));
    N(nmds) = size(MvIm,1);
    
    clear MvIm
end
close all

% Save source data files
VAL = catuneven(values,NaN);
T = table(VAL(:,1),VAL(:,2),VAL(:,3),VAL(:,4),VAL(:,5),VAL(:,6),...
    'VariableNames',names);
writetable(T,'Fig2d','Delimiter','\t','WriteRowNames',false)  

% Plot figure
f3 = figure; hold on;
set(f3,'Units','centimeters','position',[5 5 6.7 5.8],'defaultAxesFontSize',8);
plot([0 4],[0 0],':','color',[.7 .7 .7],'LineWidth',1.5);
beeswarmplot(catuneven(values(1:2:end),NaN),'orientation','left',...
    'color',colors(1:2:end),'markersize',6);
distributionPlot(values(1:2:end),'histOri','left','color',colors(1:2:end),...
    'widthDiv',[2 1],'showMM',6,'FaceAlpha',.4,'divFactor',2);
beeswarmplot(catuneven(values(2:2:end),NaN),'orientation','right',...
    'color',colors(2:2:end),'markersize',6);
distributionPlot(values(2:2:end),'histOri','right','color',colors(2:2:end),...
    'widthDiv',[2 2],'showMM',6,'FaceAlpha',.4,'divFactor',2);
xticklabels({'CA1','CA2/3','DG'});
ylabel(yl);
ylim(ylims);

%set(gca, 'YScale', 'log')

pos = get(f3,'Position');
set(f3,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(f3,sprintf('%s-SI%.2f.pdf',param,SI),'-dpdf','-r0')

kruskalwallisNdunns(catuneven(KWvalues,NaN),names);

% Print group-wise N and p-values to console
fprintf('Number of observations (N)\n');
for nmds = 1:length(names)
    fprintf('%s: %i\n',names{nmds},N(nmds));
end

fprintf('\nSigned ranks p-value\n');
for nmds = 1:length(names)
    fprintf('%s: %d\n',names{nmds},p(nmds));
end

clear param nmds mds yl names val Foutliers colors OutlierThr p N


%% FIGURE 3
%% A) Place maps color-coded by significance
names = {'CA1PV','CA1SOM','CA3PV','CA3SOM','DGPV','DGSOM'};
labels = {'CA1-PV','CA1-SOM','CA2/3-PV','CA2/3-SOM','DG-PV','DG-SOM'};
colors = {[0 1 1],[.5 .5 .5]};

% Make color bar figure
% colbar(1,:,:) = ones(256,3)-repmat(0:1/255:1,3,1)'.*colors{1};
% colbar(2,:,:) = ones(256,3)-repmat(0:1/255:1,3,1)'.*colors{2};
% figure; image(colbar);
% box off

ff = figure; sgtitle('DF/F');
set(ff,'Units','centimeters','defaultAxesFontSize',8,'Position',[5 5 8.5 14]);
fz = figure; sgtitle('z-scored');
set(fz,'Units','centimeters','defaultAxesFontSize',8,'Position',[20 5 8.5 14]);

for nmds = 1:length(names)
    eval(sprintf('mds = %s;',names{nmds}))
    gn = 0;
    for nds = 1:length(mds)
        for n = 1:length(mds{nds}.cells)
            gn = gn+1;
            for c = 1:length(mds{nds}.metadata.categories)
                ftr = nanmean(cat(2,...
                    mds{nds}.cells{n}.categories{c}.dFoY{:}),2);
                ztr = nanmean(cat(2,...
                    mds{nds}.cells{n}.categories{c}.dZoY{:}),2);
                
                % Smoothing - optional
                ftr = movmean(ftr,5,'omitnan');
                ztr = movmean(ztr,5,'omitnan');
                
                DF(:,gn,c) = ftr;
                DZ(:,gn,c) = ztr;                
                
                % Coloring
                ftr = ftr-min(ftr);
                ztr = ztr-min(ztr);
                ftr = ftr/max(ftr);
                ztr = ztr/max(ztr);
                
                if mds{nds}.cells{n}.spatial_P(1)<.05
                    ftr = ones(size(ftr,1),3)-ftr*colors{1};
                    ztr = ones(size(ztr,1),3)-ztr*colors{1};
                else                 
                    ftr = ones(size(ftr,1),3)-ftr*colors{2};
                    ztr = ones(size(ztr,1),3)-ztr*colors{2};
                end
                
                DFpl(gn,:,:,c) = ftr;
                DZpl(gn,:,:,c) = ztr;
            end
        end
    end
    
    for c = 1:length(mds{nds}.metadata.categories)
        [~,maxbin] = max(DF(:,:,c),[],1);
        [~,order{nmds}(:,c)] = sort(maxbin);
    end
    
    % Set NaN values to white color
    DZpl(isnan(DZpl)) = 1;
    DFpl(isnan(DFpl)) = 1;
    
    
    DFplot{nmds} = DFpl;
    DZplot{nmds} = DZpl;
    
    xv = (.1:.025:1.9) * 2; % length(m)
    yv = 1:gn;
    
    figure(ff)
    subplot(3,2,nmds)
    image(xv, yv, squeeze(DFpl(order{nmds}(:,1),:,:,1)));
    xticks([0 4]);
    xlabel('Track distance (m)');
    ylabel('Cell #');
    yticklabels('');
    unitstr = sprintf('N = %i',gn);
    text(0.3,gn-gn/11,unitstr)
    ylabel('Cell #');    
    title(labels{nmds});
    %caxis([-.1 .7]); %colorbar
    
    figure(fz)
    subplot(3,2,nmds)
    image(xv, yv, DZpl(order{nmds}(:,1),:,:,1));
    %caxis([-.2 1.8]);
    %colorbar;
    xlabel('Track distance (m)');
    yticklabels('');
    unitstr = sprintf('N = %i',gn);
    text(0.3,gn-gn/11,unitstr)
    %text(.3,gn-gn/11,labels{nmds});
    ylabel('Cell #');
    title(labels{nmds});
    clear DF DZ
end

% Write to PDF
pos = get(ff,'Position');
set(ff,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(ff,'DFoF_Placemaps_all,','-dpdf','-r0')

pos = get(fz,'Position');
set(fz,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(fz,'DZoZ_Placemaps_all,','-dpdf','-r0')

clear mds nds nmds n c r names gn fi fnan ftr P xv yv znan ztr maxbin ff fz order labels unitstr

%% B) Place tuning pie charts
names = {'CA1PV','CA1SOM','CA3PV','CA3SOM','DGPV','DGSOM'};
dslabels = {'CA1-PV','CA1-SOM','CA2/3-PV','CA2/3-SOM','DG-PV','DG-SOM'};
f1 = figure;%('Position',[200,100,120,800]);
set(f1,'Units','centimeters','Position',[4 4 1.2*2.54 6*2.54],'defaultAxesFontSize',12);
measure = 'spatial_P'; % 'spatial_P','positional_P'

for nmds = 1:length(names)
    eval(sprintf('mds = %s;',names{nmds}))
    spatial_P = metastatistics(mds,'cell',measure,'categories',1,...
        'rescategories',1,'plotting',false)';
    tuned = length(find(spatial_P<.05));
    untuned = length(find(spatial_P>.05));
    
    tbl(1:2,nmds) = cat(1,tuned, untuned);
    
    figure(f1);
    subplot(6,1,nmds);%,'Position',[0 1-nmds/6 .8 0.05])
    labels = {sprintf('%.1f%%',tuned*100/length(spatial_P)),...
        sprintf('%.1f%%',untuned*100/length(spatial_P))};
    pie(tbl(:,nmds),labels);
    ax = gca;
    ax.Colormap = [1 0 0; 0 0 1; .5 .5 .5];
    title(dslabels{nmds},'FontSize',8);
end

% Save source data files
PCtbl = tbl(1:2,:)./sum(tbl(1:2,:))*100;
T = table(PCtbl(1,:)',PCtbl(2,:)',...
    'VariableNames',{'Tuned%','Untuned%'},'RowNames',names);
writetable(T,'Fig3b','Delimiter','\t','WriteRowNames',true)  

set(findobj(f1,'type','text'),'fontsize',7)

% Save as pdf
pos = get(f1,'Position');
set(f1,'PaperPositionMode','Auto','PaperUnits','Inches',...
    'PaperSize',[pos(3), pos(4)],'Renderer','Painters')
print(f1,'PlacePiecharts.pdf','-dpdf','-r0')

[Chi2,p] = chi2tbltest(tbl);
fprintf('Chi2 = %.2f; p = %d\n',Chi2,p);
clear names labels Chi2 p pos ax nmds mds nds spatial_P tuned untuned tbl dslabels...
    SI

%% C) Spatial Info (per AUC) vs. bootstrap violin plots
names = {'CA1PV','CA1SOM','CA3PV','CA3SOM','DGPV','DGSOM'};
param = 'SIperAUC';%'MovImmRatio','speedmod','trialvar','sessionstab','spatial_coherence','spatialinfo','SIperAUC','positional_info'
rparam = sprintf([param '_r']);
ylims = [-.2 2];
yl = 'SI (normalized)';
OutlierThr = 3; % In interquartile ranges, use tukey procedure
colors = {[0 .5 0],[.5 1 .5],[0 0 1],[.5 .5 1],[1 0 0],[1 .5 .5]};
randcolor = [.8 .8 .5];
SI = 1;

clear values rvalues
for nmds = 1:length(names)
    eval(sprintf('mds = %s;',names{nmds}))
    clear val
    
    val(:,1) = metastatistics(mds,'cell',param,'categories',1,...
        'rescategories',1,'plotting',false,'SI',SI)';
    val(:,2) = metastatistics(mds,'cell',rparam,'categories',1,...
        'rescategories',1,'plotting',false)';
    
    val = tukeyOutlierRemoval(val,OutlierThr);
    
    p(nmds) = signrank(val(:,1),val(:,2));
    N(nmds) = size(val,1);
    
    values{nmds} = val(:,1);
    rvalues{nmds} = val(:,2);
end
close all

% Build 'global' array with all groups for stakeplot
valA = cat(1,values,rvalues);
valA = catuneven(reshape(valA,1,[]),NaN);

% Save source data files
sdfnames = cat(1,names,strcat(names,' rand'));
sdfnames = reshape(sdfnames,1,[]);
T = array2table(valA,'VariableNames',sdfnames);
writetable(T,'Fig3c','Delimiter','\t','WriteRowNames',false)  

f3 = figure('Units','centimeters','position',[4,4,6,4],'defaultAxesFontSize',8);
%title(param);
stakeplot(valA,'color',[.3 .3 .3],'LineAlpha',.2,'spread',.45,'jitter',0,'LineWidth',.5);
distributionPlot(values,'histOri','left','color',colors,'widthDiv',[2 1],'showMM',6,'FaceAlpha',.7);
distributionPlot(rvalues,'histOri','right','color',randcolor,'widthDiv',[2 2],'showMM',6,'FaceAlpha',.7);
xticks([])
%xticklabels({'CA1 PV','CA1 SOM','CA2/3 PV','CA2/3 SOM','DG PV','DG SOM'});
%xtickangle(45);
ylim(ylims);
ylabel(yl);
%set(gca, 'YScale', 'log')

% Create PDF output
% set(findobj(f3,'type','text'),'fontsize',7)
pos = get(f3,'Position');
set(f3,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(f3,sprintf('%s_vs_chance.pdf',param),'-dpdf','-painters')

% Plot KruskalWallisAnova across grops to console
kruskalwallisNdunns(catuneven(values,NaN),names);

% Print group-wise N and p-values to console
fprintf('Number of observations (N)\n');
for nmds = 1:length(names)
    fprintf('%s: %i\n',names{nmds},N(nmds));
end

fprintf('\nSigned ranks p-value\n');
for nmds = 1:length(names)
    fprintf('%s: %d\n',names{nmds},p(nmds));
end

clear param rparam rval rvalues nmds mds yl names val...
    colors oval f3 pos SI ylims N p OutlierThr randcolor values

%% D) Spatial coherence vs. bootstrap violin plots
names = {'CA1PV','CA1SOM','CA3PV','CA3SOM','DGPV','DGSOM'};
param = 'spatial_coherence';%'MovImmRatio','speedmod','trialvar','sessionstab','spatial_coherence','spatialinfo','SIperAUC'
rparam = sprintf([param '_r']);
ylims = [0 1];
yl = 'Spatial coherence';
OutlierThr = 3; % In interquartile ranges, use tukey procedure
colors = {[0 .5 0],[.5 1 .5],[0 0 1],[.5 .5 1],[1 0 0],[1 .5 .5]};
randcolor = [.8 .8 .5];
SI = 1;

clear values rvalues
for nmds = 1:length(names)
    eval(sprintf('mds = %s;',names{nmds}))
    clear val
    
    val(:,1) = metastatistics(mds,'cell',param,'categories',1,...
        'rescategories',1,'plotting',false,'SI',SI)';
    val(:,2) = metastatistics(mds,'cell',rparam,'categories',1,...
        'rescategories',1,'plotting',false)';
    
    val = tukeyOutlierRemoval(val,OutlierThr);
    
    p(nmds) = signrank(val(:,1),val(:,2));
    N(nmds) = size(val,1);
    
    values{nmds} = val(:,1);
    rvalues{nmds} = val(:,2);
end
close all

% Build 'global' array with all groups for stakeplot
valA = cat(1,values,rvalues);
valA = catuneven(reshape(valA,1,[]),NaN);

% Save source data files
sdfnames = cat(1,names,strcat(names,' rand'));
sdfnames = reshape(sdfnames,1,[]);
T = array2table(valA,'VariableNames',sdfnames);
writetable(T,'Fig3d','Delimiter','\t','WriteRowNames',false)  

f3 = figure('Units','centimeters','position',[4,4,6,4],'defaultAxesFontSize',8);
%title(param);
stakeplot(valA,'color',[.3 .3 .3],'LineAlpha',.2,'spread',.45,'jitter',0,'LineWidth',.5);
distributionPlot(values,'histOri','left','color',colors,'widthDiv',[2 1],'showMM',6,'FaceAlpha',.7);
distributionPlot(rvalues,'histOri','right','color',randcolor,'widthDiv',[2 2],'showMM',6,'FaceAlpha',.7);
xticks([]);
% xticklabels({'CA1 PV','CA1 SOM','CA3 PV','CA3 SOM','DG PV','DG SOM'});
% xtickangle(45);
ylim(ylims);
ylabel(yl);
%set(gca, 'YScale', 'log')

% Create PDF output
% set(findobj(f3,'type','text'),'fontsize',7)
pos = get(f3,'Position');
set(f3,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(f3,sprintf('%s_vs_chance.pdf',param),'-dpdf','-painters')

% Plot KruskalWallisAnova across grops to console
kruskalwallisNdunns(catuneven(values,NaN),names);

% Print group-wise N and p-values to console
fprintf('Number of observations (N)\n');
for nmds = 1:length(names)
    fprintf('%s: %i\n',names{nmds},N(nmds));
end

fprintf('\nSigned ranks p-value\n');
for nmds = 1:length(names)
    fprintf('%s: %d\n',names{nmds},p(nmds));
end

clear param rparam rval rvalues nmds mds yl names val...
    colors oval f3 pos SI ylims N p OutlierThr randcolor values

%% E) Within-session stability vs. bootstrap violin plots
names = {'CA1PV','CA1SOM','CA3PV','CA3SOM','DGPV','DGSOM'};
param = 'sessionstab';%'MovImmRatio','speedmod','trialvar','sessionstab','spatial_coherence','spatialinfo','SIperAUC'
rparam = sprintf([param '_r']);
ylims = [-.5 1];
yl = 'Within-session stability';
OutlierThr = 3; % In interquartile ranges, use tukey procedure
colors = {[0 .5 0],[.5 1 .5],[0 0 1],[.5 .5 1],[1 0 0],[1 .5 .5]};
randcolor = [.8 .8 .5];
SI = 1;

clear values rvalues
for nmds = 1:length(names)
    eval(sprintf('mds = %s;',names{nmds}))
    clear val
    
    val(:,1) = metastatistics(mds,'cell',param,'categories',1,...
        'rescategories',1,'plotting',false,'SI',SI)';
    val(:,2) = metastatistics(mds,'cell',rparam,'categories',1,...
        'rescategories',1,'plotting',false)';
    
    val = tukeyOutlierRemoval(val,OutlierThr);
    
    p(nmds) = signrank(val(:,1),val(:,2));
    N(nmds) = size(val,1);
    
    values{nmds} = val(:,1);
    rvalues{nmds} = val(:,2);
end
close all

% Build 'global' array with all groups for stakeplot
valA = cat(1,values,rvalues);
valA = catuneven(reshape(valA,1,[]),NaN);

% Save source data files
sdfnames = cat(1,names,strcat(names,' rand'));
sdfnames = reshape(sdfnames,1,[]);
T = array2table(valA,'VariableNames',sdfnames);
writetable(T,'Fig3e','Delimiter','\t','WriteRowNames',false)  

f3 = figure('Units','centimeters','position',[4,4,6,4.5],'defaultAxesFontSize',8);
%title(param);
stakeplot(valA,'color',[.3 .3 .3],'LineAlpha',.2,'spread',.45,'jitter',0,'LineWidth',.5);
distributionPlot(values,'histOri','left','color',colors,'widthDiv',[2 1],'showMM',6,'FaceAlpha',.7);
distributionPlot(rvalues,'histOri','right','color',randcolor,'widthDiv',[2 2],'showMM',6,'FaceAlpha',.7);
xticklabels({'CA1 PV','CA1 SOM','CA2/3 PV','CA2/3 SOM','DG PV','DG SOM'});
xtickangle(45);
ylim(ylims);
ylabel(yl);
%set(gca, 'YScale', 'log')

% Create PDF output
% set(findobj(f3,'type','text'),'fontsize',7)
pos = get(f3,'Position');
set(f3,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(f3,sprintf('%s_vs_chance.pdf',param),'-dpdf','-painters')

% Plot KruskalWallisAnova across grops to console
kruskalwallisNdunns(catuneven(values,NaN),names);

% Print group-wise N and p-values to console
fprintf('Number of observations (N)\n');
for nmds = 1:length(names)
    fprintf('%s: %i\n',names{nmds},N(nmds));
end

fprintf('\nSigned ranks p-value\n');
for nmds = 1:length(names)
    fprintf('%s: %d\n',names{nmds},p(nmds));
end

clear param rparam rval rvalues nmds mds yl names val...
    colors oval f3 pos SI ylims N p OutlierThr randcolor values


%% FIGURE 4
%% A,B + Suppl3h) Fam/Nov scatterplots and violin plots
% Take moving data only
signal = 'dFoT'; % 'dFoT' or 'zscored'
savefigs = true;
outlier_thr = 3; % Tukey procedure, 1.5 (mild), 3.0 (extreme), inf to disable
names = {'CA1PV','CA1SOM','CA3PV','CA3SOM','DGPV','DGSOM'};
labels = {'CA1-PV','CA1-SOM','CA2/3-PV','CA2/3-SOM','DG-PV','DG-SOM'};
colors = {[0 .5 0],[.5 1 .5],[0 0 1],[.5 .5 1],[1 0 0],[1 .5 .5]};
ncolors = {[.15 .65 .15],[.65 1 .65],[.25 .25 1],[.65 .65 1],[1 .25 .25],[1 .65 .65]};
sc_lim = [0 1; 0 1; 0 1; 0 1; 0 1; 0 .5];
box_lim = [0 1.1];


if strcmp(signal,'zscored')
    unitstr = '(z-scored)';
else
    unitstr = '(dF/F)';   
end

xl = ['Activity familiar ',unitstr];
yl = ['Activity novel ',unitstr];

f1 = figure('Position',[800 200 500 700]);
set(f1,'Units','Inches','defaultAxesFontSize',8,'Position',[2 2 3.8 5.6]);

for nmds = 1:length(names)
    eval(sprintf('mds = %s;',names{nmds}))
    gn = 0;
    for nds = 1:length(mds)
        for n = 1:length(mds{nds}.cells)
            gn = gn+1;
            for c = 1:2
                tr = cat(2,mds{nds}.cells{n}.categories{c}.(signal){:});
                mvg = cat(2,mds{nds}.metadata.categories{c}.moving{:});
                rate(gn,c) = nanmean(tr(mvg));
                %rate(gn,c) = nanmean(tr);
                
                % Introduced 221209, replace AUCrate in dataset with updated
                mds{nds}.cells{n}.AUCrate(c) = nanmean(tr(mvg));
            end
            
            mds{nds}.cells{n}.FN_ratio = rate(gn,1)/rate(gn,2);
            mds{nds}.cells{n}.ZdiffScore = abs((rate(gn,1)-rate(gn,2)));
        end
    end
    
    % retain new fields (diffScore, FN_ratio) with datasets.
    eval(sprintf('%s = mds;',names{nmds}))
    
    if ~isinf(outlier_thr)
        rate = tukeyOutlierRemoval(rate,outlier_thr);
    end
    
    Frates{nmds} = rate(:,1);
    Nrates{nmds} = rate(:,2);
    values{nmds} = (rate(:,2)-rate(:,1))';
    p(nmds) = signrank(rate(:,1),rate(:,2));
    
    subplot(3,2,nmds);
    hold on;
    plot(sc_lim(nmds,:),sc_lim(nmds,:),'--','Color',[.5 .5 .5]);
    scatter(rate(:,1),rate(:,2),5,colors{nmds},'filled');
    title(names{nmds}); %xlim(sc_lim); ylim(sc_lim);
    xlabel(xl); ylabel(yl);
    xlim(sc_lim(nmds,:)); ylim(sc_lim(nmds,:));
    xticks([-1 -.5 0 .5 1 1.5]);
    yticks([-1 -.5 0 .5 1 1.5]);
    %legend(sprintf('p = %.4f',p(nmds)));
    %legend(sprintf('p = %d',p(nmds)));
    title(labels{nmds});
    clear rate
end

% Plot Fam vs. Nov Mean z-rate (violin plot)

% Build 'global' array with all groups for stakeplot
rateA = cat(1,Frates,Nrates);
rateA = catuneven(reshape(rateA,1,[]),NaN);

% Save source data files
sdfnames = cat(1,strcat(names,' fam'),strcat(names,' nov'));
sdfnames = reshape(sdfnames,1,[]);
T = array2table(rateA,'VariableNames',sdfnames);
writetable(T,'Fig4ab','Delimiter','\t','WriteRowNames',false)

f2 = figure;%('position',[200,200,350,300]);
set(f2,'Units','Inches','defaultAxesFontSize',8,'Position',[1,1,3.2,2.6]);
stakeplot(rateA,'color',[.3 .3 .3],'LineAlpha',.2,'spread',.45,'jitter',0,'LineWidth',.5);
distributionPlot(Frates,'histOri','left','color',colors,'widthDiv',[2 1],'showMM',6,'FaceAlpha',.7);
distributionPlot(Nrates,'histOri','right','color',ncolors,'widthDiv',[2 2],'showMM',6,'FaceAlpha',.7);
ylim(box_lim);
xticklabels(labels);
xtickangle(45);
ylabel(['Mean activity rate ',unitstr]);

% Plot activity ratios (Suppl Fig 3H)
f3 = figure;
set(f3,'Units','Inches','defaultAxesFontSize',8,'Position',[1,1,3.2,2.6]);
plot([0 7],[0 0],':','color',[.7 .7 .7],'LineWidth',1.5);
beeswarmplot(catuneven(values,NaN),'color',colors,'markersize',8);
distributionPlot(values,'color',colors,'showMM',6,'FaceAlpha',.3);
xticklabels(labels);
xtickangle(45);
xlim([0 7]);
ylabel('Activity difference (nov/fam)');

% Save source data files
T = array2table(catuneven(values,NaN),'VariableNames',names);
writetable(T,'SupplFig3h','Delimiter','\t','WriteRowNames',false)

% Save figures
if savefigs
    pos = get(f1,'Position');
    set(f1,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
    print(f1,sprintf('FamNovScatter.pdf'),'-dpdf','-painters')
    
    pos = get(f2,'Position');
    set(f2,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
    print(f2,sprintf('FamNovViolins.pdf'),'-dpdf','-painters')

    pos = get(f3,'Position');
    set(f3,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
    print(f3,sprintf('FamNovDifferenceViolins.pdf'),'-dpdf','-painters')
end

% Print group-wise N and p-values to console
fprintf('Number of observations (N)\n');
for nmds = 1:length(names)
    fprintf('%s: %i\n',names{nmds},length(find(~isnan(values{nmds}))));
end

fprintf('\nSigned ranks p-value\n');
for nmds = 1:length(names)
    fprintf('%s: %d\n',names{nmds},p(nmds));
end

% KKW statistics for activity differences (Supplementary 3H)
kruskalwallisNdunns(catuneven(values,NaN),labels);

clear gn names xl yl tr mvg colors labels mds nmds nds p c pos ncolors colors Nrates Frates n...
    sc_lim SI box_lim Foutliers f1 f2 f3 foutl notOutl noutl outlier_thr savefigs signal unitstr values...
    labels valA rateA

%% C,D) SOM Axon hilus vs. ML activity plot
outlier_thr = 3; % 3 extreme, 1.5 moderate, Inf to disable
smoothFct = 5; % For duplicate removal
corrThr = 0.7; % For duplicate removal
ylims = [-.05 .5; -.05 .5];
colors = {[1 .7 .7],[1 .85 .85]};
labels = {'Fam','Nov'};
signaltype = 'dFoT_fneu'; % 'dFoT','dFoT_fneu'
printfigs = true;

% Get Hilus axon data
it = 0;
for ds = 1:length(AxonSomHilus)
    % Remove highly correlated boutons on local copy of the dataset
    HilWoD{ds} = remove_duplicates(AxonSomHilus{ds},corrThr,...
        'categories',[1 2],'signal',signaltype,'smoothFct',smoothFct);
    
    % Iterate through unique axons, get mean activity in FAM and NOV.
    for n = 1:length(HilWoD{ds}.cells)
        it = it+1;
        for c = [1 2] % Fam and Nov only
            signal = cat(2,HilWoD{ds}.cells{n}.categories{c}.(signaltype){:});
            moving = cat(2,HilWoD{ds}.metadata.categories{c}.moving{:});
            HilusDF(it,c) = nanmean(signal(moving));
        end
    end
end

% Get ML axon data
it = 0;
for ds = 1:length(AxonSomML)
    % Remove highly correlated boutons on local copy of the dataset
    MlWoD{ds} = remove_duplicates(AxonSomML{ds},corrThr,...
        'categories',[1 2],'signal',signaltype,'smoothFct',smoothFct);
    
    % Iterate through unique axons, get mean activity in FAM and NOV.
    for n = 1:length(MlWoD{ds}.cells)
        it = it+1;
        for c = [1 2] % Fam and Nov only
            signal = cat(2,MlWoD{ds}.cells{n}.categories{c}.(signaltype){:});
            moving = cat(2,MlWoD{ds}.metadata.categories{c}.moving{:});
            MlDF(it,c) = nanmean(signal(moving));
        end
    end
end

HilusDF = tukeyOutlierRemoval(HilusDF,outlier_thr);
MlDF = tukeyOutlierRemoval(MlDF,outlier_thr);

% Save source data files
T = array2table(MlDF,'VariableNames',{'ML fam','ML nov'});
writetable(T,'Fig4c','Delimiter','\t','WriteRowNames',false)

T = array2table(HilusDF,'VariableNames',{'Hilus fam','Hilus nov'});
writetable(T,'Fig4d','Delimiter','\t','WriteRowNames',false)

% =========== STAKE / VIOLIN PLOT ============
f1 = figure;
set(f1,'Units','Inches','defaultAxesFontSize',8,'Position',[1,1,1.3,4]);
subplot(2,1,1); hold on
stakeplot(HilusDF,'color',[.3 .3 .3],'LineAlpha',.2,'spread',.45,'jitter',0,'LineWidth',.5);
distributionPlot(HilusDF(:,1),'histOri','left','color',colors{1},...
    'widthDiv',[2 1],'showMM',6,'FaceAlpha',.7);
distributionPlot(HilusDF(:,2),'histOri','right','color',colors{2},...
    'widthDiv',[2 2],'showMM',6,'FaceAlpha',.7);
%title('Hilar SOM axons');
ylabel('Activity (dF/F)');
ylim(ylims(1,:));
xticks([]);
% xticks([.8 1.2]);
% xticklabels({'fam','nov'});
% xtickangle(45);
xlim([.35 1.65])

subplot(2,1,2); hold on
stakeplot(MlDF,'color',[.3 .3 .3],'LineAlpha',.2,'spread',.45,'jitter',0,'LineWidth',.5);
distributionPlot(MlDF(:,1),'histOri','left','color',colors{1},...
    'widthDiv',[2 1],'showMM',6,'FaceAlpha',.7);
distributionPlot(MlDF(:,2),'histOri','right','color',colors{2},...
    'widthDiv',[2 2],'showMM',6,'FaceAlpha',.7);
%title('ML SOM axons');
ylim(ylims(2,:))
ylabel('Activity (dF/F)');
xticks([.8 1.2]);
xticklabels({'fam','nov'});
xtickangle(45);
xlim([.35 1.65])

% =========== SCATTER PLOT ============
f2 = figure;
set(f2,'Units','Inches','defaultAxesFontSize',8,'Position',[1,1,2,4]);
subplot(2,1,1); hold on
scatter(MlDF(:,1),MlDF(:,2),5,colors{1},'filled');
plot([0 .5],[0 .5],'--','Color',[.5 .5 .5 .5]);
xlabel('Activity familiar (dF/F)'); ylabel('Activity novel (dF/F)');
%title('ML axons');
xlim([0 .5]); ylim([0 .5]);
xticks([-1 -.5 0 .5 1 1.5]);
yticks([-1 -.5 0 .5 1 1.5]);

subplot(2,1,2); hold on
scatter(HilusDF(:,1),HilusDF(:,2),5,colors{1},'filled');
plot([0 .5],[0 .5],'--','Color',[.5 .5 .5 .5]);
xlabel('Activity familiar (dF/F)'); ylabel('Activity novel (dF/F)');
%title('Hilus axons');
xlim([0 .5]); ylim([0 .5]);
xticks([-1 -.5 0 .5 1 1.5]);
yticks([-1 -.5 0 .5 1 1.5]);

if printfigs
    pos = get(f1,'Position');
    set(f1,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
    print(f1,sprintf('MLvsHilus_violins.pdf'),'-dpdf','-painters')
    
    pos = get(f2,'Position');
    set(f2,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
    print(f2,sprintf('MLvsHilus_scatters.pdf'),'-dpdf','-painters')
end

% STATISTICS
fprintf('SOM axon activitiy hilus:\n');
pairedSampleTest(HilusDF(:,1),HilusDF(:,2),'labels',labels,'verbose',true);
fprintf('SOM axon activitiy molecular layer:\n');
pairedSampleTest(MlDF(:,1),MlDF(:,2),'labels',labels,'verbose',true);

clear ds allDF outlierThr MlWoD HilWoD signal moving n it labels colors...
    MlDF HilusDF ylims c corrThr f1 outlier_thr pos signaltype smoothFct


%% Figure 5
%% C) Fraction of cells expressing h4MDi (Immuno)
PVfract = 0.94; % Aurore data, hard coded here
SOMfract = 0.97; % Aurore data, hard coded here

figure; 
subplot(2,1,1);
pie([PVfract 1-PVfract]);
ax = gca;
ax.Colormap = [0 .85 .65; .5 .5 .5];

subplot(2,1,2);
pie([SOMfract 1-SOMfract]);
ax = gca;
ax.Colormap = [.85 0 .65; .5 .5 .5];

clear PVfract SOMfract ax

%% D) Response of h4MDi+ DG PV and SOM IN to clz
% Make sure AUC rate conversion (see IN scripts) was run, should be done
% for Sep/1/2023 and following datasets!
PVrates = metastatistics(DREADDinPV,'cell','AUCrate','categories',[1 2],...
    'rescategories',[1 2],'combination','none','plotting',false);
SOMrates = metastatistics(DREADDinSOM,'cell','AUCrate','categories',[1 2],...
    'rescategories',[1 2],'combination','none','plotting',false);

PVrates = tukeyOutlierRemoval(PVrates,3);
SOMrates = tukeyOutlierRemoval(SOMrates,3);

% Save source data files
ratesA = catuneven({PVrates(:,1),PVrates(:,2),SOMrates(:,1),SOMrates(:,2)},NaN);
T = array2table(ratesA,'VariableNames',{'PV bsl','PV clz','SOM bsl','SOM clz'});
writetable(T,'Fig5d','Delimiter','\t','WriteRowNames',false)  

f1 = figure('position',[300 300 120 180]); hold on;
stakeplot(PVrates,'x',1,'color',[.3 .3 .3],'LineAlpha',.2,'spread',.45,'jitter',0,'LineWidth',.5);
distributionPlot(PVrates(:,1),'xValues',1,'histOri','left','color',[0 .55 .45],'widthDiv',[2 1],'showMM',6,'FaceAlpha',.7);
distributionPlot(PVrates(:,2),'xValues',1,'histOri','right','color',[0 .85 .65],'widthDiv',[2 2],'showMM',6,'FaceAlpha',.7);

stakeplot(SOMrates,'x',2,'color',[.3 .3 .3],'LineAlpha',.2,'spread',.45,'jitter',0,'LineWidth',.5);
distributionPlot(SOMrates(:,1),'xValues',2,'histOri','left','color',[.65 0 .85],'widthDiv',[2 1],'showMM',6,'FaceAlpha',.7);
distributionPlot(SOMrates(:,2),'xValues',2,'histOri','right','color',[.8 .5 1],'widthDiv',[2 2],'showMM',6,'FaceAlpha',.7);

ylim([-.1 1]);
xticks([1 2]);

pos = get(f1,'Position');
set(f1,'PaperPositionMode','Auto','PaperUnits','points','PaperSize',[pos(3), pos(4)])
print(f1,sprintf('INclzResponse.pdf'),'-dpdf','-painters')

% Statistics
pairedSampleTest(PVrates(:,1),PVrates(:,2),'labels',{'PVbsl','PVclz'},'verbose',true);
pairedSampleTest(SOMrates(:,1),SOMrates(:,2),'labels',{'SOMbsl','SOMclz'},'verbose',true);

% Cohen's D for paired samples
PVcohensD = (nanmean(PVrates(:,1))-nanmean(PVrates(:,2))) / nanstd(PVrates(:,1)-PVrates(:,2))
SOMcohensD = (nanmean(SOMrates(:,1))-nanmean(SOMrates(:,2))) / nanstd(SOMrates(:,1)-SOMrates(:,2))

clear PVrates SOMrates f1


%% FIGURE 6
%% A,E) Individual examples of place fields
plot_placeactivity(DREADDpv{10},106,0);
plot_placeactivity(DREADDsom{10},88,0);

%% B,F) DREADD Fraction of active and place cells (Including 2-way ANOVA)
cats = [2 4]; % Categories to compare (matters only for the stats section)
labels = {'fam bsl','nov bsl','fam clz','nov clz'};
PltOrd = [1 3 2 4]; % Order of categories in bars
colors = [0 0.447 0.741 ; 0.301 0.745 0.933 ; 0.85 0.325 0.098 ; 0.929 0.694 0.125];

% GET PV DATA
PVplcN = metastatistics(DREADDpv,'dataset','ncells','categories',1:4,...
    'rescategories',1:4,'active',1/60,'SI',.05,'plotting',false);
PVnact = metastatistics(DREADDpv,'dataset','ncells','categories',1:4,...
    'rescategories',1:4,'active',1/60,'SI',1,'plotting',false);
PVntot = metastatistics(DREADDpv,'dataset','ncells','categories',1:4,...
    'rescategories',1:4,'plotting',false);

PVplcF = PVplcN ./ PVntot .* 100; % CAVE: 'f' stands for 'fraction'
PVfact = PVnact ./ PVntot .* 100;

PVcontMtrx = [sum(PVplcN(:,cats(1))) sum(PVntot(:,cats(1)))-sum(PVplcN(:,cats(1)));...
    sum(PVplcN(:,cats(2))) sum(PVntot(:,cats(2)))-sum(PVplcN(:,cats(2)))];

% Save source data files
vnames = cat(2,strcat(labels,' % active'),strcat(labels,' % place'));
T = array2table(cat(2,PVfact,PVplcF),'VariableNames',vnames);
writetable(T,'Fig6b','Delimiter','\t','WriteRowNames',false)   

% GET SOM DATA
SOMnplc = metastatistics(DREADDsom,'dataset','ncells','categories',1:4,...
    'rescategories',1:4,'active',1/60,'SI',.05,'plotting',false);
SOMnact = metastatistics(DREADDsom,'dataset','ncells','categories',1:4,...
    'rescategories',1:4,'active',1/60,'SI',1,'plotting',false);
SOMntot = metastatistics(DREADDsom,'dataset','ncells','categories',1:4,...
    'rescategories',1:4,'plotting',false);

SOMplcF = SOMnplc ./ SOMntot .* 100;
SOMfact = SOMnact ./ SOMntot .* 100;

SOMcontMtrx = [sum(SOMnplc(:,cats(1))) sum(SOMntot(:,cats(1)))-sum(SOMnplc(:,cats(1)));...
    sum(SOMnplc(:,cats(2))) sum(SOMntot(:,cats(2)))-sum(SOMnplc(:,cats(2)))];

% Save source data files
vnames = cat(2,strcat(labels,' % active'),strcat(labels,' % place'));
T = array2table(cat(2,SOMfact,SOMplcF),'VariableNames',vnames);
writetable(T,'Fig6f','Delimiter','\t','WriteRowNames',false)   

% PLOTTING
f1 = figure;
set(f1,'Units','centimeters','defaultAxesFontSize',8,'Position',[4,4,8,4]);

subplot(1,2,1); hold on
bA=bar(sum(PVnact(:,PltOrd))./sum(PVntot).*100,'FaceColor','flat'); % CAVE: This is DIFFERENT from the average of the fractions!
bP=bar(sum(PVplcN(:,PltOrd))./sum(PVntot).*100,'FaceColor','flat');
bA.CData = colors(PltOrd,:);
bP.CData = colors(PltOrd,:);
bP.LineStyle='--';bP.LineWidth=.5; bA.EdgeColor= [.5 .5 .5];
bA.LineWidth=.5; alpha(.6); ylabel('% of cells'); xticks(1:4); xtickangle(45);

plot([1 2],PVplcF(:,PltOrd([1 2]))','-','color',[.2 .2 .2],'LineWidth',.5);
plot([3 4],PVplcF(:,PltOrd([3 4]))','-','color',[.2 .2 .2],'LineWidth',.5);

%title('DREADD PV');
set(gca,'xticklabel',labels(PltOrd),'Fontsize',8,'fontname','arial','Fontweight','normal',...
    'TickDir','out','lineWidth',.5,'box','off');

subplot(1,2,2); hold on
bA=bar(sum(SOMnact(:,PltOrd))./sum(SOMntot).*100,'FaceColor','flat');
bP=bar(sum(SOMnplc(:,PltOrd))./sum(SOMntot).*100,'FaceColor','flat');
bA.CData = colors(PltOrd,:);
bP.CData = colors(PltOrd,:);
bP.LineStyle='--';bP.LineWidth=.5; bA.EdgeColor= [.5 .5 .5];
bA.LineWidth=.5; alpha(.6); ylabel('% of cells'); xticks(1:4); xtickangle(45);

plot([1 2],SOMplcF(:,PltOrd([1 2]))','-','color',[.2 .2 .2],'LineWidth',.5);
plot([3 4],SOMplcF(:,PltOrd([3 4]))','-','color',[.2 .2 .2],'LineWidth',.5);

%title('DREADD SOM');
set(gca,'xticklabel',labels(PltOrd),'Fontsize',8,'fontname','arial','Fontweight','normal',...
    'TickDir','out','lineWidth',.5,'box','off');

% SAVE FIGURE
pos = get(f1,'Position');
set(f1,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[pos(3), pos(4)])
print(f1,sprintf('CellFractions.pdf'),'-dpdf','-painters')

% STATS
fprintf('Testing %s vs. %s\n',labels{cats(1)},labels{cats(2)})
[~,pf] = fishertest(PVcontMtrx);
[~,ptp] = ttest(PVplcF(:,cats(1)),PVplcF(:,cats(2)));
[psp,~] = signrank(PVplcF(:,cats(1)),PVplcF(:,cats(2)));
[~,pta] = ttest(PVfact(:,cats(1)),PVfact(:,cats(2)));
[psa,~] = signrank(PVfact(:,cats(1)),PVfact(:,cats(2)));
fprintf('PV Fisher p = %.6f\n',pf)
fprintf('PV Placefield fraction ttest p = %.6f signrank p = %.6f\n',ptp, psp)
fprintf('PV Active fraction ttest p = %.6f signrank p = %.6f\n',pta, psa)
fprintf('PV N sessions %i\n',size(PVplcF,1))

[~,pf] = fishertest(SOMcontMtrx);
[~,ptp] = ttest(SOMplcF(:,cats(1)),SOMplcF(:,cats(2)));
[psp,~] = signrank(SOMplcF(:,cats(1)),SOMplcF(:,cats(2)));
[~,pta] = ttest(SOMfact(:,cats(1)),SOMfact(:,cats(2)));
[psa,~] = signrank(SOMfact(:,cats(1)),SOMfact(:,cats(2)));
fprintf('SOM Fisher p = %.6f\n',pf)
fprintf('SOM Placefield fraction ttest p = %.6f signrank p = %.6f\n',ptp, psp)
fprintf('SOM Active fraction ttest p = %.6f signrank p = %.6f\n',pta, psa)
fprintf('SOM N sessions %i\n',size(SOMplcF,1))

% Between genotypes comparison
[~,pFam] = ttest2(PVplcF(:,1),SOMplcF(:,1));
[~,pNov] = ttest2(PVplcF(:,2),SOMplcF(:,2));
fprintf('\nComparing genotypes SOM vs. PV (t-test, bsl): Fam p = %.4f, Nov p = %.4f\n',pFam,pNov);
fprintf('Comparing genotypes SOM vs. PV (t-test, bsl): Fam p = %.4f, Nov p = %.4f\n',...
    ranksum(PVplcF(:,1),SOMplcF(:,1)), ranksum(PVplcF(:,2),SOMplcF(:,2)));

% 2-way ANOVA (PV)
FNlabels = repmat([0 1 0 1],size(PVplcF,1),1);
BClabels = repmat([0 0 1 1],size(PVplcF,1),1);
[~,tbl,stats] = anovan(PVplcF(:),{FNlabels(:),BClabels(:)},...
    "Model","interaction","Varnames",["IsNov","IsClz"]);
[results,~,~,gnames] = multcompare(stats,"Dimension",[1 2]);

fprintf('2-way ANOVA with interactions for PV-h4MDi, Fam-Nov, Bsl-Clz:\n');
fprintf('%s\t%s\t%s\t%s\t%s\t%s\t%s\n',tbl{1,1},tbl{1,2},tbl{1,3},tbl{1,4},tbl{1,5},tbl{1,6},tbl{1,7})
for n = 2:size(tbl,1)
    fprintf('%s\t%.2e\t%i\t%.2g\t%.2f\t%.2g\t%.4g\n',tbl{n,1},tbl{n,2},tbl{n,3},tbl{n,4},tbl{n,5},tbl{n,6},tbl{n,7})
end

% 2-way ANOVA (SOM)
FNlabels = repmat([0 1 0 1],size(SOMplcF,1),1);
BClabels = repmat([0 0 1 1],size(SOMplcF,1),1);
[~,tbl,stats] = anovan(SOMplcF(:),{FNlabels(:),BClabels(:)},...
    "Model","interaction","Varnames",["IsNov","IsClz"]);
[results,~,~,gnames] = multcompare(stats,"Dimension",[1 2]);

fprintf('2-way ANOVA with interactions for SOM-h4MDi, Fam-Nov, Bsl-Clz:\n');
fprintf('%s\t%s\t%s\t%s\t%s\t%s\t%s\n',tbl{1,1},tbl{1,2},tbl{1,3},tbl{1,4},tbl{1,5},tbl{1,6},tbl{1,7})
for n = 2:size(tbl,1)
    fprintf('%s\t%.2e\t%i\t%.2g\t%.2f\t%.2g\t%.4g\n',tbl{n,1},tbl{n,2},tbl{n,3},tbl{n,4},tbl{n,5},tbl{n,6},tbl{n,7})
end

clear PVnplc PVntot PVfplc SOMnplc SOMntot SOMfplc cats PVcontMtrx SOMcontMtrx...
    PVnact PVfact SOMnact SOMfact colors PltOrd labels nSOMcells nPVcells...
    pf ptp psp pta psa outlier_thr bA bP pFam pNov

%% C,D,G,H) Transientrate and SI violin plots DREADD dataset -- new w/ stakeplot
% This version is built for pair-wise comparisons primarily.

% SETTING PARAMETERS
plotcats = [1 3; 2 4]; % These are the pairs of categories querried for the violin plots, 
% dim 1 defines number of violins on each subplot
selcats = {[1 3], [2 4]}; % These are the categories used for cell selections 
% (set e.g. to [1 2] to get cells meeting criteria during baseline, [1:6]
% for cells meeting crietria under any condition).
fields = {'transientrate','spatialinfo'}; % This can be any property that is 
% a subfield of 'cells' in the Data structure, e.g 'transientrate','spatial_P',etc.
ACTthr = 1/60;
pSIthr = .05;
pPFthr = 1;
savefigs = false;
combination = 'or';
outlierthr = 3; % Threshold (in interquartile ranges) for tukey outlier removal, set to Inf to disable feature
colors = {[0 0.447 0.741],[0.301 0.745 0.933],...
    [0.85 0.325 0.098],[0.929 0.694 0.125],...
    [.5 .65 .5], [.7 .85 .7]};
xlabels = {'fam-bsl','nov-bsl','fam-clz','nov-clz'};
dst = .2; % Distance of xticklabels
ylabels = {'Transient rate (Hz)','spatial info (bits/sec)'}; % Needs to match fields
ylims = [-.01 .1; -.1 1.5]; % Plotting range for each parameter. Needs to match fields
sourcenames = {'Fig6c','Fig6d';'Fig6g','Fig6h'};


% EXTRACTING DATA
for f = 1:length(fields)
    % Filling in values for each group of categories.
    for cs = 1:size(plotcats,1)
        PVvalues = metastatistics(DREADDpv,'cell',fields{f},'categories',selcats{cs},'SI',pSIthr,...
            'PF',pPFthr,'active',ACTthr,'combination',combination,'rescategories',plotcats(cs,:),'plotting',false);
        SOMvalues = metastatistics(DREADDsom,'cell',fields{f},'categories',selcats{cs},'SI',pSIthr,...
            'PF',pPFthr,'active',ACTthr,'combination',combination,'rescategories',plotcats(cs,:),'plotting',false);
        
        % Optional: outlier removal here
        PVvalues = tukeyOutlierRemoval(PVvalues,outlierthr);
        SOMvalues = tukeyOutlierRemoval(SOMvalues,outlierthr);

        PVres.(fields{f})(cs,[1 2]) = {PVvalues(:,1),PVvalues(:,2)};
        SOMres.(fields{f})(cs,[1 2]) = {SOMvalues(:,1),SOMvalues(:,2)};
    end
end

% PLOTTING
f1 = figure;
set(f1,'Units','centimeters','defaultAxesFontSize',8,'Position',[4,4,7,8.8]);

for f = 1:length(fields) 
   % Build 'global' array with all groups for PV and SOM animals
   PVa = catuneven(reshape(PVres.(fields{f})',1,[]),NaN);
   SOMa = catuneven(reshape(SOMres.(fields{f})',1,[]),NaN);
   
   % Save source data files
   T = array2table(PVa,'VariableNames',xlabels(plotcats'));
   writetable(T,sourcenames{1,f},'Delimiter','\t','WriteRowNames',false)
   T = array2table(SOMa,'VariableNames',xlabels(plotcats'));
   writetable(T,sourcenames{2,f},'Delimiter','\t','WriteRowNames',false)
   
   % Plot PV data
   subplot(2,length(fields),f); hold on;
   stakeplot(PVa,'color',[.3 .3 .3],'LineAlpha',.2,'spread',.45,'jitter',0,'LineWidth',.5);
   distributionPlot(PVres.(fields{f})(:,1),'histOri','left','color',colors(plotcats(:,1)),'widthDiv',[2 1],'showMM',6,'FaceAlpha',.7);
   distributionPlot(PVres.(fields{f})(:,2),'histOri','right','color',colors(plotcats(:,2)),'widthDiv',[2 2],'showMM',6,'FaceAlpha',.7);
   
   title('PV Cre');
   ylim(ylims(f,:));
   ylabel(ylabels{f});
   xticks(sort([(1:size(plotcats,1))-dst (1:size(plotcats,1))+dst]));
   xticklabels(xlabels(plotcats'));
   xtickangle(45);
   xlim([.4 length(fields)+.6]);
   
   % Plot SOM data
   subplot(2,length(fields),length(fields)+f);
   stakeplot(SOMa,'color',[.3 .3 .3],'LineAlpha',.2,'spread',.45,'jitter',0,'LineWidth',.5);
   distributionPlot(SOMres.(fields{f})(:,1),'histOri','left','color',colors(plotcats(:,1)),'widthDiv',[2 1],'showMM',6,'FaceAlpha',.7);
   distributionPlot(SOMres.(fields{f})(:,2),'histOri','right','color',colors(plotcats(:,2)),'widthDiv',[2 2],'showMM',6,'FaceAlpha',.7);
   
   title('SOM Cre');
   ylim(ylims(f,:));
   ylabel(ylabels{f});
   xticks(sort([(1:size(plotcats,1))-dst (1:size(plotcats,1))+dst]));
   xticklabels(xlabels(plotcats'));
   xtickangle(45);
   xlim([.4 length(fields)+.6]);
end

% SAVE FIGURE
if savefigs
    pos = get(f1,'Position');
    set(f1,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[pos(3), pos(4)])
    print(f1,sprintf('SInRate.pdf'),'-dpdf','-painters')
end

% STATISTICS
for f = 1:length(fields)
    fprintf('PV Cells %s:\n',fields{f});
    for cs = 1:size(plotcats,1)
        pairedSampleTest(PVres.(fields{f}){cs,1},PVres.(fields{f}){cs,2},...
            'verbose',true,'labels',xlabels(plotcats(cs,:)));
    end
    
    fprintf('SOM Cells %s:\n',fields{f});
    for cs = 1:size(plotcats,1)
        pairedSampleTest(SOMres.(fields{f}){cs,1},SOMres.(fields{f}){cs,2},...
            'verbose',true,'labels',xlabels(plotcats(cs,:)));
    end
end

clear pPFthr pSIthr ACTthr cats colors dst plotcats outlierthr fields...
    ylabels ylims f cs f1 PVvalues SOMvalues pPVs pPVt pSOMs pSOMt...
    PVres SOMres xlabels selcats combination PVa SOMa

%% I) DREADD mCherry Control Fraction of active and place cells
cats = [2 4]; % Categories to compare (matters only for the stats section)
labels = {'fam bsl','nov bsl','fam clz','nov clz'};
PltOrd = [1 3 2 4]; % Order of categories in bars
colors = [0 0.447 0.741 ; 0.301 0.745 0.933 ; 0.85 0.325 0.098 ; 0.929 0.694 0.125];

% GET PV DATA
CTRplcN = metastatistics(DREADDctr,'dataset','ncells','categories',1:4,...
    'rescategories',1:4,'active',1/60,'SI',.05,'plotting',false);
CTRnact = metastatistics(DREADDctr,'dataset','ncells','categories',1:4,...
    'rescategories',1:4,'active',1/60,'SI',1,'plotting',false);
CTRntot = metastatistics(DREADDctr,'dataset','ncells','categories',1:4,...
    'rescategories',1:4,'plotting',false);

CTRplcF = CTRplcN ./ CTRntot .* 100; % CAVE: 'f' stands for 'fraction'
CTRfact = CTRnact ./ CTRntot .* 100;

CTRcontMtrx = [sum(CTRplcN(:,cats(1))) sum(CTRntot(:,cats(1)))-sum(CTRplcN(:,cats(1)));...
    sum(CTRplcN(:,cats(2))) sum(CTRntot(:,cats(2)))-sum(CTRplcN(:,cats(2)))];

% Save source data files
vnames = cat(2,strcat(labels,' % active'),strcat(labels,' % place'));
T = array2table(cat(2,CTRfact,CTRplcF),'VariableNames',vnames);
writetable(T,'Fig6i','Delimiter','\t','WriteRowNames',false)   

% PLOTTING
f1 = figure; hold on
set(f1,'Units','centimeters','defaultAxesFontSize',8,'Position',[4,4,4,4]);

bA=bar(sum(CTRnact(:,PltOrd))./sum(CTRntot).*100,'FaceColor','flat'); % CAVE: This is DIFFERENT from the average of the fractions!
bP=bar(sum(CTRplcN(:,PltOrd))./sum(CTRntot).*100,'FaceColor','flat');
bA.CData = colors(PltOrd,:);
bP.CData = colors(PltOrd,:);
bP.LineStyle='--';bP.LineWidth=.5; bA.EdgeColor= [.5 .5 .5];
bA.LineWidth=.5; alpha(.6); ylabel('% of cells'); xticks(1:4); xtickangle(45);

plot([1 2],CTRplcF(:,PltOrd([1 2]))','-','color',[.2 .2 .2],'LineWidth',.5);
plot([3 4],CTRplcF(:,PltOrd([3 4]))','-','color',[.2 .2 .2],'LineWidth',.5);

set(gca,'xticklabel',labels(PltOrd),'Fontsize',8,'fontname','arial','Fontweight','normal',...
    'TickDir','out','lineWidth',.5,'box','off');

% SAVE FIGURE
pos = get(f1,'Position');
set(f1,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[pos(3), pos(4)])
print(f1,sprintf('Ctr CellFractions.pdf'),'-dpdf','-painters')

% STATS
fprintf('Testing %s vs. %s\n',labels{cats(1)},labels{cats(2)})
[~,pf] = fishertest(CTRcontMtrx);
[~,ptp] = ttest(CTRplcF(:,cats(1)),CTRplcF(:,cats(2)));
[psp,~] = signrank(CTRplcF(:,cats(1)),CTRplcF(:,cats(2)));
[~,pta] = ttest(CTRfact(:,cats(1)),CTRfact(:,cats(2)));
[psa,~] = signrank(CTRfact(:,cats(1)),CTRfact(:,cats(2)));
fprintf('CTR Fisher p = %.6f\n',pf)
fprintf('CTR Placefield fraction ttest p = %.6f signrank p = %.6f\n',ptp, psp)
fprintf('CTR Active fraction ttest p = %.6f signrank p = %.6f\n',pta, psa)
fprintf('CTR N sessions %i\n',size(CTRplcF,1))

clear CTRnplc CTRntot CTRfplc cats CTRcontMtrx...
    CTRnact CTRfact colors PltOrd labels nCTRcells...
    pf ptp psp pta psa outlier_thr bA bP pFam pNov

%% J,K) Rate and spatialinfo Control mCherry injection BSL vs CLZ
% This version is built for pair-wise comparisons primarily.

% SETTING PARAMETERS
plotcats = [1 3; 2 4]; % These are the pairs of categories querried for the violin plots, 
% dim 1 defines number of violins on each subplot
selcats = {[1 3], [2 4]}; % These are the categories used for cell selections 
% (set e.g. to [1 2] to get cells meeting criteria during baseline, [1:6]
% for cells meeting crietria under any condition).
fields = {'transientrate','spatialinfo'}; % This can be any property that is 
% a subfield of 'cells' in the Data structure, e.g 'transientrate','spatial_P',etc.
ACTthr = 1/60;
pSIthr = .05;
pPFthr = 1;
combination = 'or';
outlierthr = 3; % Threshold (in interquartile ranges) for tukey outlier removal, set to Inf to disable feature
colors = {[0 0.447 0.741],[0.301 0.745 0.933],...
    [0.85 0.325 0.098],[0.929 0.694 0.125],...
    [.5 .65 .5], [.7 .85 .7]};
xlabels = {'fam-bsl','nov-bsl','fam-clz','nov-clz'};
dst = .2; % Distance of xticklabels
ylabels = {'Transient rate (Hz)','spatial info (bits/sec)'}; % Needs to match fields
ylims = [-.01 .15; -.08 1.5]; % Plotting range for each parameter. Needs to match fields
sourcenames = {'Fig6j','Fig6k'};


% EXTRACTING DATA
for f = 1:length(fields)
    % Filling in values for each group of categories.
    for cs = 1:size(plotcats,1)
        SOMvalues = metastatistics(DREADDctr,'cell',fields{f},'categories',selcats{cs},'SI',pSIthr,...
            'PF',pPFthr,'active',ACTthr,'combination',combination,'rescategories',plotcats(cs,:),'plotting',false);
        
        % Optional: outlier removal here
        SOMvalues = tukeyOutlierRemoval(SOMvalues,outlierthr);

        SOMres.(fields{f})(cs,[1 2]) = {SOMvalues(:,1),SOMvalues(:,2)};
    end
end

% PLOTTING
f1 = figure;
%set(f1,'Units','centimeters','defaultAxesFontSize',8,'Position',[4,8,7,6]);
set(f1,'Units','centimeters','defaultAxesFontSize',8,'Position',[4,4,7,4.4]);

for f = 1:length(fields) 
   % Build 'global' array with all groups for PV and SOM animals
   SOMa = catuneven(reshape(SOMres.(fields{f})',1,[]),NaN);
   
   % Save source data files
   T = array2table(SOMa,'VariableNames',xlabels(plotcats'));
   writetable(T,sourcenames{f},'Delimiter','\t','WriteRowNames',false)
   
   subplot(1,length(fields),f);
   stakeplot(SOMa,'color',[.2 .2 .2],'LineAlpha',.1,'spread',.45,'jitter',0,'LineWidth',.5);
   distributionPlot(SOMres.(fields{f})(:,1),'histOri','left','color',colors(plotcats(:,1)),'widthDiv',[2 1],'showMM',6,'FaceAlpha',.7);
   distributionPlot(SOMres.(fields{f})(:,2),'histOri','right','color',colors(plotcats(:,2)),'widthDiv',[2 2],'showMM',6,'FaceAlpha',.7);
   %title('SOM Cre');
   ylim(ylims(f,:));
   ylabel(ylabels{f});
   xticks(sort([(1:size(plotcats,1))-dst (1:size(plotcats,1))+dst]));
   xticklabels(xlabels(plotcats'));
   xtickangle(45);
   xlim([.4 length(fields)+.6]);
end

% SAVE FIGURE
pos = get(f1,'Position');
set(f1,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[pos(3), pos(4)])
print(f1,sprintf('SInRate_ctr.pdf'),'-dpdf','-painters')

% STATISTICS
for f = 1:length(fields)
    fprintf('Clozapine control %s:\n',fields{f});
    for cs = 1:size(plotcats,1)
        pairedSampleTest(SOMres.(fields{f}){cs,1},SOMres.(fields{f}){cs,2},...
            'verbose',true,'labels',xlabels(plotcats(cs,:)));
    end
end

clear pPFthr pSIthr ACTthr cats colors dst plotcats outlierthr fields...
    ylabels ylims f cs f1 PVvalues SOMvalues pPVs pPVt pSOMs pSOMt...
    PVres SOMres xlabels selcats combination PVa SOMa

%% L,M) 1-Way Kruskal-Wallis for PV, SOM, Control + 2-way ANOVA
% SETTING PARAMETERS
plotcats = [1 3; 2 4]; % These are the pairs of categories querried for the violin plots, 
% dim 1 defines number of violins on each subplot
selcats = {[1 3], [2 4]}; % These are the categories used for cell selections 
% (set e.g. to [1 2] to get cells meeting criteria during baseline, [1:6]
% for cells meeting crietria under any condition).
fields = {'transientrate','spatialinfo'}; % This can be any subfield of 'cells', e.g 'transientrate','spatial_P',etc.
datasets = {'DREADDpv','DREADDsom','DREADDctr'};
ACTthr = 1/60;
pSIthr = .05;
pPFthr = 1;
savefigs = false;
combination = 'or';
outlierthr = 3; % Threshold (in interquartile ranges) for tukey outlier removal, set to Inf to disable feature
colors = repmat({[0.85 0.325 0.098],[0.929 0.694 0.125]},3,1);
sourcenames = {'Fig6l','Fig6m'};

res = {[],[]};

% EXTRACTING DATA
for f = 1:length(fields)
    for ds = 1:length(datasets)
        eval(sprintf('metads = %s;',datasets{ds}));
        
        for cs = 1:size(plotcats,1)
            values = metastatistics(metads,'cell',fields{f},'categories',selcats{cs},'SI',pSIthr,...
                'PF',pPFthr,'active',ACTthr,'combination',combination,'rescategories',plotcats(cs,:),'plotting',false);
            
            ratios = values(:,2)./values(:,1);
            %ratios = real(log2(ratios));
            
            % Optional: outlier removal here
            ratios = tukeyOutlierRemoval(ratios,outlierthr);
            
            resCell{f,cs,ds} = ratios;
            
            % ANOVA labels (column 2 = PV/SOM/Ctr, 3 = Fam/Nov, 4 =
            % compound
            ratios(:,2) = ds-1; % Column 2: PV = 0, SOM = 1, Ctr = 2
            ratios(:,3) = cs-1; % Column 3: 0 = fam, 1 = nov
            ratios(:,4) = ratios(:,2) + 3*(cs-1); % 012 = fam, 345 = nov
            
            res{f} = cat(1,res{f},ratios);
        end
    end
    
    % PLOT RATIOS
    f1 = figure; hold on;
    plotvals = real(log2(res{f}(:,1)));
    [~,isoutlier] = tukeyOutlierRemoval(plotvals,3);
    
    set(f1,'Units','centimeters','defaultAxesFontSize',10,'Position',[4,4,8.6,6.6]);
    
    % Group index = PVfam(0), PVnov(1), SOMfam(0+2), SOMnov(1+2);
    distributionPlot(plotvals(~isoutlier),'groups',...
        res{f}(~isoutlier,4),...
        'color',colors,'showMM',6,'distWidth',.9,...
        'addSpread',0,'histOpt',1,'divFactor',3,'FaceAlpha',.5);
    
    plotSpread(plotvals(~isoutlier),'distributionIdx',...
        res{f}(~isoutlier,4),...
        'spreadFcn',{'lin',10},'spreadWidth',.9,...
        'distributionColors',[0 0 0],...
        'distributionMarkerSize',1,'binWidth',.05);
        
    ylabel(sprintf('Log2 %s ratio clz/bsl',fields{f}));
    xticks(0:5);
    xticklabels({'PV Fam','SOM Fam','Ctr Fam','PV Nov','SOM Nov','Ctr Nov'});
    xtickangle(45);
    
    pos = get(f1,'Position');
    set(f1,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[pos(3), pos(4)])
    print(f1,sprintf('%s ratios 3way.pdf',fields{f}),'-dpdf','-painters')    
    
    % Kruskal-Wallis
    fprintf('%s familiar kruskal wallis and Dunn''s',fields{f})
    kruskalwallisNdunns(catuneven(resCell(f,1,:),NaN),{'PV','SOM','Ctr'});
    fprintf('%s novel kruskal wallis and Dunn''s',fields{f})
    kruskalwallisNdunns(catuneven(resCell(f,2,:),NaN),{'PV','SOM','Ctr'});
    
    % Save source data files
    T = array2table(catuneven(cat(3,resCell(f,1,:),resCell(f,2,:)),NaN),'VariableNames',...
        {'PV fam','SOM fam','Ctr fam','PV nov','SOM nov','Ctr nov'});
    writetable(T,sourcenames{f},'Delimiter','\t','WriteRowNames',false)
    
    % 2-WAY ANOVA
    % Limit to experimental groups here
    isexp = res{f}(:,2)<2;
    [~,tbl,stats] = anovan(res{f}(isexp,1),{res{f}(isexp,2),res{f}(isexp,3)},...
        "Model","interaction","Varnames",["Strain","Nov"]);
    fprintf('2-way ANOVA with interactions for %s in SOM-PV, FAM-NOV:\n',fields{f});
    fprintf('%s\t%s\t%s\t%s\t%s\t%s\t%s\n',tbl{1,1},tbl{1,2},tbl{1,3},tbl{1,4},tbl{1,5},tbl{1,6},tbl{1,7})
    for n = 2:size(tbl,1)
        fprintf('%s\t%.2e\t%i\t%.2g\t%.2f\t%.2g\t%.4g\n',tbl{n,1},tbl{n,2},tbl{n,3},tbl{n,4},tbl{n,5},tbl{n,6},tbl{n,7})
    end
    
    figure;
    [compc,~,~,gnames] = multcompare(stats,"Dimension",[1 2]);
    fprintf('\nDunn\''s post-hoc test\n')
    fprintf('%s\t%s\t%s\n','Group1','Group2','p')
    for n = 1:size(compc,1)
        fprintf('%s\t%s\t%.2e\n',gnames{compc(n,1)},gnames{compc(n,2)},compc(n,6))
    end
end

clear pPFthr pSIthr ACTthr cats colors dst plotcats outlierthr fields...
    ylabels ylims f cs f1 values ratios...
    xlabels selcats combination PVa SOMa metads

%% Text only: F-statistic to assess for unequal variances (pre-post clozapine)
% F = ratio mean variance clz / mean variance bsl - larger variance goes to
% numerator!; Degrees of freedom = sample size (N) -1, respectively. Take
% p*2 for 2-tailed test!
%
% p = (1-fpdf(F,N1-1,N2-1))*2;

cats = [1 3]; % [2 4] for novel bsl vs. clz

PVplcF = metastatistics(DREADDpv,'cell','transientrate','categories',cats,'SI',.05,...
    'PF',1,'active',1/60,'combination','or','rescategories',cats,'plotting',false);
SOMplcF = metastatistics(DREADDsom,'cell','transientrate','categories',cats,'SI',.05,...
    'PF',1,'active',1/60,'combination','or','rescategories',cats,'plotting',false);

FpvF = var(PVplcF(:,2))/var(PVplcF(:,1)); if FpvF<1; FpvF = 1/FpvF; end
FsomF = var(SOMplcF(:,2))/var(SOMplcF(:,1)); if FsomF<1; FsomF = 1/FsomF; end

pPV = (1-fcdf(FpvF,size(PVplcF,1)-1,size(PVplcF,1)-1))*2;
pSOM = (1-fcdf(FsomF,size(SOMplcF,1)-1,size(SOMplcF,1)-1))*2;

fprintf('Comparison categories %i vs. %i\n',cats(1),cats(2));
fprintf('PV: F = %.2f, p = %.2e, %i degrees of freedom\n',FpvF,pPV,size(PVplcF,1)-1)
fprintf('SOM: F = %.2f, p = %.2e, %i degrees of freedom\n',FsomF,pSOM,size(SOMplcF,1)-1)

clear PVfplc SOMfplc cats


%% FIGURE 7
%% A,C,E) Place maps
plot_placemaps_meta(DREADDpv,1,'SI',.05,'active',1/60,'clims',[0 .2],'figscale',[.5 .5]);
input('Press any key to continue'); close all;
plot_placemaps_meta(DREADDpv,3,'SI',.05,'active',1/60,'clims',[0 .2],'figscale',[.5 .5]);
input('Press any key to continue'); close all;
plot_placemaps_meta(DREADDsom,1,'SI',.05,'active',1/60,'clims',[0 .2],'figscale',[.5 .5]);
input('Press any key to continue'); close all;
plot_placemaps_meta(DREADDsom,3,'SI',.05,'active',1/60,'clims',[0 .2],'figscale',[.5 .5]);
input('Press any key to continue'); close all;
plot_placemaps_meta(DREADDctr,1,'categories',[1 2],'SI',.05,'active',1/60,...
    'combination','or','clims',[0 .2],'figscale',[.5 .5]);
input('Press any key to continue'); close all;
plot_placemaps_meta(DREADDctr,3,'categories',[3 4],'SI',.05,'active',1/60,...
    'combination','or','clims',[0 .2],'figscale',[.5 .5]);
input('Press any key to continue'); close all;

clear c dataset_name2 dFoYplot sortindices axis_value

%% B) Place-field correlations PV pre- and post clozapine violin plots
labels = {'bsl FN','clz FN'};
colors = {[0 0.447 0.741],[0.85 0.325 0.098]};

[res,~] = metacorr_pairwise(DREADDpv,'categories',1:4,'active',1/60,...
    'SI',.05,'combination','or','plotresults',false);

FNbsl = squeeze(res(1,2,:));
FNclz = squeeze(res(3,4,:));

f1 = figure; hold on;
set(f1,'Units','centimeters','defaultAxesFontSize',6,'Position',[4,4,2.8,4.5]);

%figure('position',[200 200 150 220]); hold on;
distributionPlot({FNbsl,FNclz},'color',colors,'showMM',6,'distWidth',.9,...
    'addSpread',0,'histOpt',1,'divFactor',3,'FaceAlpha',.5);
scatter(randn(length(FNbsl),1)/10+1,FNbsl,8,max(colors{1}-.2,0),'.')
scatter(randn(length(FNclz),1)/10+2,FNclz,8,max(colors{2}-.2,0),'.')
ylim([-.7 1.2]);
ylabel('Spatial Correlation (R)','FontSize',6);
xticklabels(labels);
xtickangle(45);
xlim([.5 2.5])
set(gca,'FontSize',6);

% Save source data files
T = array2table(catuneven({FNbsl,FNclz},NaN),'VariableNames',...
    {'FN corr bsl','FN corr clz'});
writetable(T,'Fig7b','Delimiter','\t','WriteRowNames',false)

% SAVE FIGURE
pos = get(f1,'Position');
set(f1,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[pos(3), pos(4)])
print(f1,sprintf('PV_FN_bslVSclz.pdf'),'-dpdf','-painters','-bestfit')

fprintf('Comparing PV cells with place fields during bsl, or clz respectively:\n')
unpairedSampleTest(FNbsl,FNclz,'verbose',true,'labels',labels);

clear res labels FNbsl FNclz

%% D) Place-field correlations SOM pre- and post clozapine violin plots
labels = {'bsl FN','clz FN'};
colors = {[0 0.447 0.741],[0.85 0.325 0.098]};

[res,~] = metacorr_pairwise(DREADDsom,'categories',1:4,'active',1/60,...
    'SI',.05,'combination','or','plotresults',false);

FNbsl = squeeze(res(1,2,:));
FNclz = squeeze(res(3,4,:));

f1 = figure; hold on;
set(f1,'Units','centimeters','defaultAxesFontSize',6,'Position',[4,4,2.8,4.5]);

distributionPlot({FNbsl,FNclz},'color',colors,'showMM',6,'distWidth',.9,...
    'addSpread',0,'histOpt',1,'divFactor',3,'FaceAlpha',.5);
scatter(randn(length(FNbsl),1)/10+1,FNbsl,8,max(colors{1}-.2,0),'.')
scatter(randn(length(FNclz),1)/10+2,FNclz,8,max(colors{2}-.2,0),'.')
ylim([-.7 1.2]);
ylabel('Spatial Correlation (R)');
xticklabels(labels);
xtickangle(45);
set(gca,'FontSize',6);

% Save source data files
T = array2table(catuneven({FNbsl,FNclz},NaN),'VariableNames',...
    {'FN corr bsl','FN corr clz'});
writetable(T,'Fig7d','Delimiter','\t','WriteRowNames',false)

% SAVE FIGURE
pos = get(f1,'Position');
set(f1,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[pos(3), pos(4)])
print(f1,sprintf('SOM_FN_bslVSclz.pdf'),'-dpdf','-painters','-bestfit')

fprintf('Comparing SOM cells with place fields during bsl, or clz respectively:\n')
unpairedSampleTest(FNbsl,FNclz,'verbose',true,'labels',labels);

clear res labels FNbsl FNclz

%% F) Place-field correlations Ctr pre- and post clozapine violin plots
labels = {'bsl FN','clz FN'};
colors = {[0 0.447 0.741],[0.85 0.325 0.098]};

[res,~] = metacorr_pairwise(DREADDctr,'categories',1:4,'active',1/60,...
    'SI',.05,'combination','or','plotresults',false);

FNbsl = squeeze(res(1,2,:));
FNclz = squeeze(res(3,4,:));

f1 = figure; hold on;
set(f1,'Units','centimeters','defaultAxesFontSize',6,'Position',[4,4,2.4,5]);

%figure('position',[200 200 150 220]); hold on;
distributionPlot({FNbsl,FNclz},'color',colors,'showMM',6,'distWidth',.9,...
    'addSpread',0,'histOpt',1,'divFactor',3,'FaceAlpha',.5);
scatter(randn(length(FNbsl),1)/10+1,FNbsl,8,max(colors{1}-.2,0),'.')
scatter(randn(length(FNclz),1)/10+2,FNclz,8,max(colors{2}-.2,0),'.')
ylim([-.7 1.2]);
ylabel('Spatial Correlation (R)','FontSize',6);
xticklabels(labels);
xtickangle(45);
xlim([.5 2.5])
set(gca,'FontSize',6);

% Save source data files
T = array2table(catuneven({FNbsl,FNclz},NaN),'VariableNames',...
    {'FN corr bsl','FN corr clz'});
writetable(T,'Fig7f','Delimiter','\t','WriteRowNames',false)

% SAVE FIGURE
pos = get(f1,'Position');
set(f1,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[pos(3), pos(4)])
print(f1,sprintf('Ctr_FN_bslVSclz.pdf'),'-dpdf','-painters','-bestfit')

fprintf('Comparing granule cells with place fields during bsl, or clz respectively:\n')
unpairedSampleTest(FNbsl,FNclz,'verbose',true,'labels',labels);

clear res labels FNbsl FNclz


%% SUPPLEMENTARY FIGURE 1 -- CELL EXAMPLES
%% D-I) Representative examples All Regions
names = {'CA1PV','CA1SOM','CA3PV','CA3SOM','DGPV','DGSOM'};
labels = {'CA1-PV','CA1-SOM','CA2/3-PV','CA2/3-SOM','DG-PV','DG-SOM'};
ex = [1 12 32; 2 7 12; 3 4 4; 4 2 1; 5 2 2; 6 7 15];
signal = {'dFoT','dFodY','dFoY'}; % {'zscored','dZodY','dZoY'}
colmap = jet;
clims = [0 .3; 0 .7; 0 .9; 0 2.5; 0 1; 0 .3];
xlims = [500 700; 0 200; 250 450; 1200 1400; 850 1050; 900 1100];
%colmap(1,:) = [1 1 1];
% SOMslow = [6 7 15; 
% PVfast = [1 1 1; 1 1 3; 1 6 1;
%col = {'b','c'};
f1 = figure;
set(f1,'defaultAxesFontSize',8,'Units','Inches','Position',[2 2 8.3 3.9]);
VarSaveStrings = strcat('SupplFig1',{'d','e','f','g','h','i'});

for e = 1:size(ex,1)
    eval(sprintf('mds = %s;',names{ex(e,1)}))
    
    % SUBPLOT 1: TRACES
    %subplot(size(ex,1),4,[(e-1)*4+1,((e-1)*4+2)])
    %subplot(size(ex,1),3,3*(e-1)+1);
    subplot(3,size(ex,1),e);
    title(labels{ex(e,1)});
    hold on;
    fr = mds{ex(e,2)}.metadata.categories{1}.acquisition_rate(1);
    
    % Put together all traces
    [tr,v,y,dZoY] = deal([]);
    famI = 0;
    
    for r = 1:length(mds{ex(e,2)}.metadata.categories{1}.y)
        famL = length(mds{ex(e,2)}.cells{ex(e,3)}.categories{1}.(signal{1}){r});
        dtr = cat(2,mds{ex(e,2)}.cells{ex(e,3)}.categories{1}.(signal{1}){r},...
            mds{ex(e,2)}.cells{ex(e,3)}.categories{2}.(signal{1}){r});
        tr = cat(2,tr,dtr);
        v = cat(2,v,mds{ex(e,2)}.metadata.categories{1}.v{r},...
            mds{ex(e,2)}.metadata.categories{2}.v{r});
        y = cat(2,y,mds{ex(e,2)}.metadata.categories{1}.y{r},...
            mds{ex(e,2)}.metadata.categories{2}.y{r});
        
        % Plot rectangle for new track
        rectangle('Position',[famI/fr,-10,famL/fr,15],'FaceColor',[.7 .7 .7],...
            'LineWidth',1,'EdgeColor',[.7 .7 .7])
        
        famI = famI + length(dtr);
    end
    clear famI
    
    t = (0:length(v)-1)/fr;
    v(v>70) = 70; % Speed-scale: max = 30 cm/s
    tr = movmean(tr,10);
    
    % Calcium trace (black)
    plot(t,1.5*tr/clims(e,2),'k');
    % Speed (blue)
    plot(t,v/20-3,'b');
    % Track location (y; red)
    plot(t,y-5,'r');
    
    %xlim([0 180]);
    xlim(xlims(e,:));
    xlabel('Time (s)');
    ylim([-5 5])
    
    % SUBPLOT 2 SPEED FIT
    %subplot(size(ex,1),3,3*(e-1)+2);
    subplot(3,size(ex,1),size(ex,1)+e);
    hold on
    xv = 0:15;
    pro = nanmean(cat(2,mds{ex(e,2)}.cells{ex(e,3)}.categories{1}.(signal{2}){:}),2);
    prosem = nanstd(cat(2,mds{ex(e,2)}.cells{ex(e,3)}.categories{1}.(signal{2}){:}),[],2)./...
        sqrt(length(mds{ex(e,2)}.cells{ex(e,3)}.categories{1}.(signal{2})));
    sR = mds{ex(e,2)}.cells{ex(e,3)}.speed_R(1);
    sM = mds{ex(e,2)}.cells{ex(e,3)}.speedmod(1);
    sp = mds{ex(e,2)}.cells{ex(e,3)}.speed_P(1);
    errorbar(xv,pro(1:16),prosem(1:16),'k');
    P = polyfit(xv(~isnan(pro(1:16))),pro(~isnan(pro(1:16)))',1);
    title(sprintf('M=%.2d, p=%.2d',sM,sp),'FontSize',6);
    plot(xv,xv*P(1)+P(2),'.','color','b');
    ylabel('Activity (DF/F)');
    xlabel('Speed (cm/s)');
    
    % Save source data files
    T = table(xv',pro(1:16),prosem(1:16),'VariableNames',{'Speed (cm/s)','Mean Activity','SEM Activity'});
    writetable(T,VarSaveStrings{e},'Delimiter','\t');
    
    % SUBPLOT 3: PLACE MAP
    %subplot(size(ex,1),3,3*(e-1)+3);
    subplot(3,size(ex,1),2*size(ex,1)+e);
    % Variable still named dZoY, but can also be DF if set appropriately
    dZoY = cat(2,dZoY,mds{ex(e,2)}.cells{ex(e,3)}.categories{1}.(signal{3}){:}); 
    dZoY = movmean(dZoY,3,1);
    xv = 0:4;
    yv = 1:size(dZoY,2);
    
    imagesc(xv,yv,dZoY');
    caxis(clims(e,:));
    colormap(colmap)
    set(gca,'box','off');
    xlabel('Track distance (m)');
    ylabel('Trial #');
    
end

pos = get(f1,'Position');
set(f1,'PaperPositionMode','Auto','PaperUnits','Inches',...
    'PaperSize',[pos(3), pos(4)],'Renderer','Painters')
print(f1,'ExampleCellsSup','-dpdf','-r0')

clear names v t tr pro P xv ax mds e ex prosem famI famL fr dtr r y yv...
    dZoY f1 labels pos sM sp sR colmap


%% SUPPLEMENTARY FIGURE 2 -- RESULTS PER ANIMAL AND DATASETS
%% Preparation: Get universal colormap across datasets
rng(0);
names = {'CA1PV','CA1SOM','CA3PV','CA3SOM','DGPV','DGSOM','DREADDpv','DREADDsom'};
animals = [];
for nmds = 1:length(names)
    eval(sprintf('mds = %s;',names{nmds}))
    for nds = 1:length(mds)
        animals(end+1) = mds{nds}.animal;
    end
end
lut(:,1) = unique(animals);
lut(:,2:4) = rand(size(lut,1),3);

clear animals mds names nmds nds

%% A) Speed modulation slope (2C)
clear animalRes sessionRes
names = {'CA1PV','CA1SOM','CA3PV','CA3SOM','DGPV','DGSOM'};
labels = {'CA1-PV','CA1-SOM','CA3-PV','CA3-SOM','DG-PV','DG-SOM'};
parameter = 'speedmod';

f1 = figure('position',[400 400 300 240]); hold on;
for nmds = 1:length(names)
    eval(sprintf('mds = %s;',names{nmds}))
    
    % Get cell-by-cell data metrics with session and animal info
    [data,cellinfo] = metastatistics(mds,'cell',parameter,...
        'categories',1,'rescategories',1,'plotting',false);
    
    % Remove extreme outliers, in line with other analyses
    [data,isoutlier] = tukeyOutlierRemoval(data,3);
    cellinfo = cellinfo(1,~isoutlier);
    
    % Plot by-animal and by-session data, get results in structures
    [session_data, animal_data] = plot_by_animal(...
        data, cellinfo, nmds, lut);
    
    animalRes{nmds} = [animal_data(:).data];
    sessionRes{nmds} = [session_data(:).data];
end

% Markup
xticks(1:6);
xticklabels(labels);
xtickangle(45);
ylabel('Speed modulation slope');

% SAVE FIGURE
pos = get(f1,'Position');
set(f1,'PaperPositionMode','Auto','PaperUnits','points','PaperSize',[pos(3), pos(4)])
print(f1,sprintf([parameter '.pdf']),'-dpdf','-painters')

% Concatenate means for session and animals
animalResA = catuneven(animalRes,NaN);
sessionResA = catuneven(sessionRes,NaN);

% Save source data files
sdfnames = cat(1,strcat(names,' animal'),strcat(names,' session'));
sdfnames = reshape(sdfnames',1,[]);
T = array2table(catuneven([animalRes,sessionRes],NaN),'VariableNames',sdfnames);
writetable(T,'SupplFig2a','Delimiter','\t','WriteRowNames',false)

% Statistics
fprintf('Animal comparison (Kruskal Wallis)\n');
kruskalwallisNdunns(animalResA,names);

fprintf('Session comparison (Kruskal Wallis)\n');
kruskalwallisNdunns(sessionResA,names);

clear mds nmds names data cellinfo animal_data session_data...
    animalRes sessionRes parameter f1 labels pos

%% B) Activity ratio moving - immobile (2D)
clear animalRes sessionRes
names = {'CA1PV','CA1SOM','CA3PV','CA3SOM','DGPV','DGSOM'};
labels = {'CA1-PV','CA1-SOM','CA3-PV','CA3-SOM','DG-PV','DG-SOM'};
parameter = 'MovImmRatio';
yaxistitle = 'Activity Ratio moving/immobile';

f1 = figure('position',[400 400 300 240]); hold on;
for nmds = 1:length(names)
    eval(sprintf('mds = %s;',names{nmds}))
    
    % Get cell-by-cell data metrics with session and animal info
    [data,cellinfo] = metastatistics(mds,'cell',parameter,...
        'categories',1,'rescategories',1,'plotting',false);
    
%     % Log2 conversion
%     data = real(log2(data));
    
    % Remove extreme outliers, in line with other analyses
    [data,isoutlier] = tukeyOutlierRemoval(data,3);
    cellinfo = cellinfo(1,~isoutlier);
    
    % Plot by-animal and by-session data, get results in structures
    [session_data, animal_data] = plot_by_animal(...
        data, cellinfo, nmds, lut);
    
    animalRes{nmds} = [animal_data(:).data];
    sessionRes{nmds} = [session_data(:).data];
end

% Markup
xticks(1:6);
xticklabels(labels);
xtickangle(45);
ylabel(yaxistitle);

% SAVE FIGURE
pos = get(f1,'Position');
set(f1,'PaperPositionMode','Auto','PaperUnits','points','PaperSize',[pos(3), pos(4)])
print(f1,sprintf([parameter '.pdf']),'-dpdf','-painters')

% Concatenate means for session and animals
animalResA = catuneven(animalRes,NaN);
sessionResA = catuneven(sessionRes,NaN);

% Save source data files
sdfnames = cat(1,strcat(names,' animal'),strcat(names,' session'));
sdfnames = reshape(sdfnames',1,[]);
T = array2table(catuneven([animalRes,sessionRes],NaN),'VariableNames',sdfnames);
writetable(T,'SupplFig2b','Delimiter','\t','WriteRowNames',false)

% Statistics
fprintf('Animal comparison (Kruskal Wallis)\n');
kruskalwallisNdunns(animalResA,names);

fprintf('Session comparison (Kruskal Wallis)\n');
kruskalwallisNdunns(sessionResA,names);

clear mds nmds names data cellinfo animal_data session_data...
    animalRes sessionRes parameter labels f1 pos

%% C) Spatial Info per AUC vs. chance (3C)
% This will be 6 individual plots (row) scatter vs. chance
names = {'CA1PV','CA1SOM','CA3PV','CA3SOM','DGPV','DGSOM'};
labels = {'CA1-PV','CA1-SOM','CA2/3-PV','CA2/3-SOM','DG-PV','DG-SOM'};
parameter = 'SIperAUC';
axistitle = 'Spatial Information';
lims = [0 2];
allSourceData = {};
SourceDataNames = {};

f1 = figure('position',[100 100 1280 184]);
for nmds = 1:length(names)
    subplot(1,length(names),nmds); hold on;
    eval(sprintf('mds = %s;',names{nmds}))
    
    % Get cell-by-cell data metrics with session and animal info
    [data,cellinfo] = metastatistics(mds,'cell',parameter,...
        'categories',1,'rescategories',1,'plotting',false);
    [rdata,~] = metastatistics(mds,'cell',[parameter '_r'],...
        'categories',1,'rescategories',1,'plotting',false);
    cdata = cat(1,rdata',data');
    
    % Remove extreme outliers, in line with other analyses
    [cdata,isoutlier] = tukeyOutlierRemoval(cdata',3); cdata = cdata';
    cellinfo = cellinfo(1,~isoutlier);
    
    % Plot by-animal and by-session data, get results in structures
    [session_data, animal_data] = scatter_by_animal(...
        cdata, cellinfo, lut);
    
    ylabel(axistitle);
    xlabel([axistitle ' (rand)']);
    title(labels{nmds});
    
    xlim(lims); ylim(lims);
    plot(lims,lims,'-.','color',[.5 .5 .5]);
    vrange = lims(2)-lims(1);
    
    % Statistics
    fprintf('%s (animals):\n',names{nmds});
    data = [animal_data(:).data];
    pA(nmds) = pairedSampleTest(data(1,:),data(2,:),'labels',{axistitle,'bootstrap'},'verbose',true);
    t1 = text(lims(1)+0.05*vrange, lims(2)-0.05*vrange,...
        sprintf('p=%.4f, N=%i',pA(nmds),size(data,2)));
    t1.FontSize = 8;
    
    printf('%s (sessions): ',names{nmds});
    data = [session_data(:).data];
    pS(nmds) = pairedSampleTest(data(1,:),data(2,:),'labels',{'bootstrap',axistitle},'verbose',true);
    t2 = text(lims(1)+0.4*vrange, lims(1)+0.075*vrange,...
        sprintf('p=%.4f, N=%i',pS(nmds),size(data,2)));
    t2.FontSize = 8;
    
    % Make array for source data
    % animal data
    SourceDataNames{end+1} = strcat([names{nmds} ' animal'],{' rand',' data',' rand SEM',' data SEM',' N cells'});
    allSourceData{end+1} = catuneven({[animal_data(:).data]',[animal_data(:).SEM]',[animal_data(:).n]},NaN);
    % session data
    SourceDataNames{end+1} = strcat([names{nmds} ' session'],{' rand',' data',' rand SEM',' data SEM',' N cells'});
    allSourceData{end+1} = catuneven({[session_data(:).data]',[session_data(:).SEM]',[session_data(:).n]},NaN);

end

% Save source data -- CAVE, order is bootstrap(1) - realData (2).
T = array2table(catuneven(allSourceData,NaN),'VariableNames',cat(2,SourceDataNames{:}));
writetable(T,'SupplFig2c','Delimiter','\t','WriteRowNames',false)

% Print stats to console
fprintf('%s per animal:\n',axistitle);
for nmds = 1:length(labels)
    fprintf('%s data vs. shuffle: p = %.4g\n',labels{nmds},pA(nmds));
end
fprintf('%s per session:\n',axistitle);
for nmds = 1:length(labels)
    fprintf('%s data vs. shuffle: p = %.4g\n',labels{nmds},pS(nmds));
end

% SAVE FIGURE
pos = get(f1,'Position');
set(f1,'PaperPositionMode','Auto','PaperUnits','points','PaperSize',[pos(3), pos(4)])
print(f1,sprintf([parameter '.pdf']),'-dpdf','-painters')

clear mds nmds names data cellinfo animal_data session_data cdata...
    animalRes sessionRes parameter axistitle labels pA pS rdata vrange...
    t1 t2 f1 isoutlier lims pos axistitle

%% D) Spatial coherence vs. chance (3D)
% This will be 6 individual plots (row) scatter vs. chance
names = {'CA1PV','CA1SOM','CA3PV','CA3SOM','DGPV','DGSOM'};
labels = {'CA1-PV','CA1-SOM','CA2/3-PV','CA2/3-SOM','DG-PV','DG-SOM'};
parameter = 'spatial_coherence';
axistitle = 'Spatial Coherence';
lims = [0.2 1];
allSourceData = {};
SourceDataNames = {};

f1 = figure('position',[100 100 1280 184]);
for nmds = 1:length(names)
    subplot(1,length(names),nmds); hold on;
    eval(sprintf('mds = %s;',names{nmds}))
    
    % Get cell-by-cell data metrics with session and animal info
    [data,cellinfo] = metastatistics(mds,'cell',parameter,...
        'categories',1,'rescategories',1,'plotting',false);
    [rdata,~] = metastatistics(mds,'cell',[parameter '_r'],...
        'categories',1,'rescategories',1,'plotting',false);
    cdata = cat(1,rdata',data');
    
    % Remove extreme outliers, in line with other analyses
    [cdata,isoutlier] = tukeyOutlierRemoval(cdata',3); cdata = cdata';
    cellinfo = cellinfo(1,~isoutlier);
    
    % Plot by-animal and by-session data, get results in structures
    [session_data, animal_data] = scatter_by_animal(...
        cdata, cellinfo, lut);
    
    ylabel(axistitle);
    xlabel([axistitle ' (rand)']);
    title(labels{nmds});
    
    xlim(lims); ylim(lims);
    plot(lims,lims,'-.','color',[.5 .5 .5]);
    vrange = lims(2)-lims(1);
    
    % Statistics
    fprintf('%s (animals):\n',names{nmds});
    data = [animal_data(:).data];
    pA(nmds) = pairedSampleTest(data(1,:),data(2,:),'labels',{axistitle,'bootstrap'},'verbose',true);
    t1 = text(lims(1)+0.05*vrange, lims(2)-0.05*vrange,...
        sprintf('p=%.4f, N=%i',pA(nmds),size(data,2)));
    t1.FontSize = 8;
    
    printf('%s (sessions): ',names{nmds});
    data = [session_data(:).data];
    pS(nmds) = pairedSampleTest(data(1,:),data(2,:),'labels',{axistitle,'bootstrap'},'verbose',true);
    t2 = text(lims(1)+0.4*vrange, lims(1)+0.075*vrange,...
        sprintf('p=%.4f, N=%i',pS(nmds),size(data,2)));
    t2.FontSize = 8;
    
    % Make array for source data
    % animal data
    SourceDataNames{end+1} = strcat([names{nmds} ' animal'],{' rand',' data',' rand SEM',' data SEM',' N cells'});
    allSourceData{end+1} = catuneven({[animal_data(:).data]',[animal_data(:).SEM]',[animal_data(:).n]},NaN);
    % session data
    SourceDataNames{end+1} = strcat([names{nmds} ' session'],{' rand',' data',' rand SEM',' data SEM',' N cells'});
    allSourceData{end+1} = catuneven({[session_data(:).data]',[session_data(:).SEM]',[session_data(:).n]},NaN);

end

% Save source data -- CAVE, order is bootstrap(1) - realData (2).
T = array2table(catuneven(allSourceData,NaN),'VariableNames',cat(2,SourceDataNames{:}));
writetable(T,'SupplFig2d','Delimiter','\t','WriteRowNames',false)
clear SourceDataNames allSourceData

% Print stats to console
fprintf('%s per animal:\n',axistitle);
for nmds = 1:length(labels)
    fprintf('%s data vs. shuffle: p = %.4g\n',labels{nmds},pA(nmds));
end
fprintf('%s per session:\n',axistitle);
for nmds = 1:length(labels)
    fprintf('%s data vs. shuffle: p = %.4g\n',labels{nmds},pS(nmds));
end

% SAVE FIGURE
pos = get(f1,'Position');
set(f1,'PaperPositionMode','Auto','PaperUnits','points','PaperSize',[pos(3), pos(4)])
print(f1,sprintf([parameter '.pdf']),'-dpdf','-painters')

clear mds nmds names data cellinfo animal_data session_data cdata...
    animalRes sessionRes parameter axistitle labels pA pS rdata vrange...
    t1 t2 f1 isoutlier lims pos axistitle

%% E) Sessionstab vs. chance (3E)
% This will be 6 individual plots (row) scatter vs. chance
names = {'CA1PV','CA1SOM','CA3PV','CA3SOM','DGPV','DGSOM'};
labels = {'CA1-PV','CA1-SOM','CA2/3-PV','CA2/3-SOM','DG-PV','DG-SOM'};
parameter = 'sessionstab';
axistitle = 'Within-session stability';
lims = [-.5 .75];
allSourceData = {};
SourceDataNames = {};

f1 = figure('position',[100 100 1280 184]);
for nmds = 1:length(names)
    subplot(1,length(names),nmds); hold on;
    eval(sprintf('mds = %s;',names{nmds}))
    
    % Get cell-by-cell data metrics with session and animal info
    [data,cellinfo] = metastatistics(mds,'cell',parameter,...
        'categories',1,'rescategories',1,'plotting',false);
    [rdata,~] = metastatistics(mds,'cell',[parameter '_r'],...
        'categories',1,'rescategories',1,'plotting',false);
    cdata = cat(1,rdata',data');
    
    % Remove extreme outliers, in line with other analyses
    [cdata,isoutlier] = tukeyOutlierRemoval(cdata',3); cdata = cdata';
    cellinfo = cellinfo(1,~isoutlier);
    
    % Plot by-animal and by-session data, get results in structures
    [session_data, animal_data] = scatter_by_animal(...
        cdata, cellinfo, lut);
    
    ylabel(axistitle);
    xlabel([axistitle ' (rand)']);
    title(labels{nmds});
    
    xlim(lims); ylim(lims);
    plot(lims,lims,'-.','color',[.5 .5 .5]);
    vrange = lims(2)-lims(1);
    
    % Statistics
    fprintf('%s (animals):\n',names{nmds});
    data = [animal_data(:).data];
    pA(nmds) = pairedSampleTest(data(1,:),data(2,:),'labels',{axistitle,'bootstrap'},'verbose',true);
    t1 = text(lims(1)+0.05*vrange, lims(2)-0.05*vrange,...
        sprintf('p=%.4f, N=%i',pA(nmds),size(data,2)));
    t1.FontSize = 8;
    
    printf('%s (sessions): ',names{nmds});
    data = [session_data(:).data];
    pS(nmds) = pairedSampleTest(data(1,:),data(2,:),'labels',{axistitle,'bootstrap'},'verbose',true);
    t2 = text(lims(1)+0.4*vrange, lims(1)+0.075*vrange,...
        sprintf('p=%.4f, N=%i',pS(nmds),size(data,2)));
    t2.FontSize = 8;
    
    % Make array for source data
    % animal data
    SourceDataNames{end+1} = strcat([names{nmds} ' animal'],{' rand',' data',' rand SEM',' data SEM',' N cells'});
    allSourceData{end+1} = catuneven({[animal_data(:).data]',[animal_data(:).SEM]',[animal_data(:).n]},NaN);
    % session data
    SourceDataNames{end+1} = strcat([names{nmds} ' session'],{' rand',' data',' rand SEM',' data SEM',' N cells'});
    allSourceData{end+1} = catuneven({[session_data(:).data]',[session_data(:).SEM]',[session_data(:).n]},NaN);

end

% Save source data -- CAVE, order is bootstrap(1) - realData (2).
T = array2table(catuneven(allSourceData,NaN),'VariableNames',cat(2,SourceDataNames{:}));
writetable(T,'SupplFig2e','Delimiter','\t','WriteRowNames',false)
clear SourceDataNames allSourceData

% Print stats to console
fprintf('%s per animal:\n',axistitle);
for nmds = 1:length(labels)
    fprintf('%s data vs. shuffle: p = %.4g\n',labels{nmds},pA(nmds));
end
fprintf('%s per session:\n',axistitle);
for nmds = 1:length(labels)
    fprintf('%s data vs. shuffle: p = %.4g\n',labels{nmds},pS(nmds));
end

% SAVE FIGURE
pos = get(f1,'Position');
set(f1,'PaperPositionMode','Auto','PaperUnits','points','PaperSize',[pos(3), pos(4)])
print(f1,sprintf([parameter '.pdf']),'-dpdf','-painters')

clear mds nmds names data cellinfo animal_data session_data cdata...
    animalRes sessionRes parameter axistitle labels pA pS rdata vrange...
    t1 t2 f1 isoutlier lims pos axistitle

%% F) Activity fam vs. nov (4B)
% Scatter fam vs. nov -- RUN 4B CODE FIRST FOR UPDATED AUCrates!
names = {'CA1PV','CA1SOM','CA3PV','CA3SOM','DGPV','DGSOM'};
labels = {'CA1-PV','CA1-SOM','CA2/3-PV','CA2/3-SOM','DG-PV','DG-SOM'};
lims = [.1 .6; .1 .6; 0 1; 0 .8; 0 .8; .1 .4];
allSourceData = {};
SourceDataNames = {};

f1 = figure('position',[100 100 1280 184]);
for nmds = 1:length(names)
    subplot(1,length(names),nmds); hold on;
    eval(sprintf('mds = %s;',names{nmds}))
    
    % Get cell-by-cell data metrics with session and animal info
    [data,cellinfo] = metastatistics(mds,'cell','AUCrate',...
        'categories',[1 2],'rescategories',[1 2],'plotting',false);
    
    % Remove extreme outliers, in line with other analyses
    [data,isoutlier] = tukeyOutlierRemoval(data,3);
    cellinfo = cellinfo(1,~isoutlier);
    
    % Plot by-animal and by-session data, get results in structures
    [session_data, animal_data] = scatter_by_animal(...
        data', cellinfo, lut);
    
    ylabel('Activity novel (dF/F)');
    xlabel('Activity familiar (dF/F)');
    title(labels{nmds});
    
    xlim(lims(nmds,:)); ylim(lims(nmds,:));
    plot(lims(nmds,:),lims(nmds,:),'-.','color',[.5 .5 .5]);
    vrange = lims(nmds,2)-lims(nmds,1);
    
    % Statistics
    fprintf('%s (animals):\n',names{nmds});
    data = [animal_data(:).data];
    pA(nmds) = pairedSampleTest(data(1,:),data(2,:),'labels',{'familiar','novel'},'verbose',true);
    t1 = text(lims(nmds,1)+0.05*vrange, lims(nmds,2)-0.05*vrange,...
        sprintf('p=%.4f, N=%i',pA(nmds),size(data,2)));
    t1.FontSize = 8;
    
    printf('%s (sessions): ',names{nmds});
    data = [session_data(:).data];
    pS(nmds) = pairedSampleTest(data(1,:),data(2,:),'labels',{'familiar','novel'},'verbose',true);
    t2 = text(lims(nmds,1)+0.4*vrange, lims(nmds,1)+0.075*vrange,...
        sprintf('p=%.4f, N=%i',pS(nmds),size(data,2)));
    t2.FontSize = 8;
    
    % Make array for source data
    % animal data
    SourceDataNames{end+1} = strcat([names{nmds} ' animal'],{' fam',' nov',' fam SEM',' nov SEM',' N cells'});
    allSourceData{end+1} = catuneven({[animal_data(:).data]',[animal_data(:).SEM]',[animal_data(:).n]},NaN);
    % session data
    SourceDataNames{end+1} = strcat([names{nmds} ' session'],{' fam',' nov',' fam SEM',' nov SEM',' N cells'});
    allSourceData{end+1} = catuneven({[session_data(:).data]',[session_data(:).SEM]',[session_data(:).n]},NaN);

end

% Save source data -- CAVE, order is bootstrap(1) - realData (2).
T = array2table(catuneven(allSourceData,NaN),'VariableNames',cat(2,SourceDataNames{:}));
writetable(T,'SupplFig2f','Delimiter','\t','WriteRowNames',false)
clear SourceDataNames allSourceData

% Print stats to console
fprintf('Fam vs. nov per animal:\n');
for nmds = 1:length(labels)
    fprintf('%s data vs. shuffle: p = %.4g\n',labels{nmds},pA(nmds));
end
fprintf('Fam vs. nov per session:\n');
for nmds = 1:length(labels)
    fprintf('%s data vs. shuffle: p = %.4g\n',labels{nmds},pS(nmds));
end

% SAVE FIGURE
pos = get(f1,'Position');
set(f1,'PaperPositionMode','Auto','PaperUnits','points','PaperSize',[pos(3), pos(4)])
print(f1,'Fam vs Nov.pdf','-dpdf','-painters')

clear mds nmds names data cellinfo animal_data session_data cdata...
    animalRes sessionRes parameter axistitle labels pA pS rdata vrange...
    t1 t2 f1 isoutlier lims pos axistitle

%% G-J) Transient rate and SI, PV:DREADD (5 D,E)
% SETTING PARAMETERS
labels = {'Transient Rate fam','Transient Rate nov','Spatial Info fam','Spatial Info nov'};
plotcats = [1 3; 2 4; 1 3; 2 4]; % These are the pairs of categories querried for the violin plots, 
% dim 1 defines number of violins on each subplot
selcats = {[1 3], [2 4], [1 3], [2 4]}; % These are the categories used for cell selections 
% (set e.g. to [1 2] to get cells meeting criteria during baseline, [1:6]
% for cells meeting crietria under any condition).
fields = {'transientrate','transientrate','spatialinfo','spatialinfo'}; % This can be any property that is 
% a subfield of 'cells' in the Data structure, e.g 'transientrate','spatial_P',etc.
ACTthr = 1/60;
pSIthr = .05;
pPFthr = 1;
combination = 'or';
lims = [0 .05; 0 .05; -.1 1; -.1 1]; % Plotting range for each parameter. Needs to match fields
sdfnames = strcat('SupplFig2',{'g','h','i','j'});

f1 = figure('position',[100 100 1000 250]); sgtitle('PV-Cre::DREADD');

for f = 1:length(labels)
    subplot(1,4,f); hold on;
    % Scatter saline vs. expt.
    [data,cellinfo] = metastatistics(DREADDpv,'cell',fields{f},'categories',selcats{f},'SI',pSIthr,...
        'PF',pPFthr,'active',ACTthr,'combination',combination,'rescategories',plotcats(f,:),'plotting',false);

    % Remove extreme outliers, in line with other analyses
    [data,isoutlier] = tukeyOutlierRemoval(data,3);
    cellinfo = cellinfo(1,~isoutlier);
    
    % Plot by-animal and by-session data, get results in structures
    [session_data, animal_data] = scatter_by_animal(...
        data', cellinfo, lut);
    
    ylabel('clozapine');
    xlabel('baseline');
    title(labels{f});
    
    xlim(lims(f,:)); ylim(lims(f,:));
    plot(lims(f,:),lims(f,:),'-.','color',[.5 .5 .5]);
    vrange = lims(f,2)-lims(f,1);
    
    % Statistics
    fprintf('PV-Cre DREADD %s (animals):\n',labels{f});
    data = [animal_data(:).data];
    pA(f) = pairedSampleTest(data(1,:),data(2,:),'labels',{'baseline','clozapine'},'verbose',true);
    t1 = text(lims(f,1)+0.05*vrange, lims(f,2)-0.05*vrange,...
        sprintf('p=%.4f, N=%i',pA(f),size(data,2)));
    %t1.FontSize = 8;
    
    printf('PV-Cre DREADD %s (sessions): ',labels{f});
    data = [session_data(:).data];
    pS(f) = pairedSampleTest(data(1,:),data(2,:),'labels',{'baseline','clozapine'},'verbose',true);
    t2 = text(lims(f,1)+0.4*vrange, lims(f,1)+0.075*vrange,...
        sprintf('p=%.4f, N=%i',pS(f),size(data,2)));
    %t2.FontSize = 8;

    % Save source data
    sdfvarnames{1} = strcat('PV-h4MDi animal',{' bsl',' clz',' bsl SEM',' clz SEM',' N cells'});
    sdfvarnames{2} = strcat('PV-h4MDi session',{' bsl',' clz',' bsl SEM',' clz SEM',' N cells'});
    T = array2table(catuneven({[animal_data(:).data]',[animal_data(:).SEM]',[animal_data(:).n],...
        [session_data(:).data]',[session_data(:).SEM]',[session_data(:).n]},NaN),...
        'VariableNames',cat(2,sdfvarnames{:}));
    writetable(T,sdfnames{f},'Delimiter','\t','WriteRowNames',false)
end

% Print stats to console
fprintf('DREADD effects per animal:\n');
for nmds = 1:length(labels)
    fprintf('%s data vs. shuffle: p = %.4g\n',labels{nmds},pA(nmds));
end
fprintf('DREADD effects per session:\n');
for nmds = 1:length(labels)
    fprintf('%s data vs. shuffle: p = %.4g\n',labels{nmds},pS(nmds));
end

% SAVE FIGURE
pos = get(f1,'Position');
set(f1,'PaperPositionMode','Auto','PaperUnits','points','PaperSize',[pos(3), pos(4)])
print(f1,sprintf(['PV-DREADD' '.pdf']),'-dpdf','-painters')

clear mds nmds names data cellinfo animal_data session_data...
    animalRes sessionRes parameter axistitle labels pA pS rdata vrange...
    ACTthr combination f fields isoutlier lims plotcats pPFthr pSIthr restable...
    selcats ylabels

%% L-O) Transient rate and SI, SOM:DREADD (5 G,H)
% SETTING PARAMETERS
labels = {'Transient Rate (Hz) fam','Transient Rate (Hz) nov','Spatial Info (bits/sec) fam','Spatial Info (bits/sec) nov'};
plotcats = [1 3; 2 4; 1 3; 2 4]; % These are the pairs of categories querried for the violin plots, 
% dim 1 defines number of violins on each subplot
selcats = {[1 3], [2 4], [1 3], [2 4]}; % These are the categories used for cell selections 
% (set e.g. to [1 2] to get cells meeting criteria during baseline, [1:6]
% for cells meeting crietria under any condition).
fields = {'transientrate','transientrate','spatialinfo','spatialinfo'}; % This can be any property that is 
% a subfield of 'cells' in the Data structure, e.g 'transientrate','spatial_P',etc.
ACTthr = 1/60;
pSIthr = .05;
pPFthr = 1;
combination = 'or';
lims = [0 .08; 0 .08; -.1 .8; -.1 .8]; % Plotting range for each parameter. Needs to match fields
sdfnames = strcat('SupplFig2',{'l','m','n','o'});

f1 = figure('position',[100 100 1000 250]); sgtitle('SOM-Cre::DREADD');

for f = 1:length(labels)
    subplot(1,4,f); hold on;
    % Scatter saline vs. expt.
    [data,cellinfo] = metastatistics(DREADDsom,'cell',fields{f},'categories',selcats{f},'SI',pSIthr,...
        'PF',pPFthr,'active',ACTthr,'combination',combination,'rescategories',plotcats(f,:),'plotting',false);

    % Remove extreme outliers, in line with other analyses
    [data,isoutlier] = tukeyOutlierRemoval(data,3);
    cellinfo = cellinfo(1,~isoutlier);
    
    % Plot by-animal and by-session data, get results in structures
    [session_data, animal_data] = scatter_by_animal(...
        data', cellinfo, lut);
    
    ylabel('clozapine');
    xlabel('baseline');
    title(labels{f});
    
    xlim(lims(f,:)); ylim(lims(f,:));
    plot(lims(f,:),lims(f,:),'-.','color',[.5 .5 .5]);
    vrange = lims(f,2)-lims(f,1);
    
    % Statistics
    fprintf('SOM-Cre DREADD %s (animals):\n',labels{f});
    data = [animal_data(:).data];
    pA(f) = pairedSampleTest(data(1,:),data(2,:),'labels',{'baseline','clozapine'},'verbose',true);
    text(lims(f,1)+0.05*vrange, lims(f,2)-0.05*vrange,...
        sprintf('p=%.4f, N=%i',pA(f),size(data,2)));
    
    printf('SOM-Cre DREADD %s (sessions): ',labels{f});
    data = [session_data(:).data];
    pS(f) = pairedSampleTest(data(1,:),data(2,:),'labels',{'baseline','clozapine'},'verbose',true);
    text(lims(f,1)+0.4*vrange, lims(f,1)+0.075*vrange,...
        sprintf('p=%.4f, N=%i',pS(f),size(data,2)));
    
    % Save source data
    sdfvarnames{1} = strcat('SOM-h4MDi animal',{' bsl',' clz',' bsl SEM',' clz SEM',' N cells'});
    sdfvarnames{2} = strcat('SOM-h4MDi session',{' bsl',' clz',' bsl SEM',' clz SEM',' N cells'});
    T = array2table(catuneven({[animal_data(:).data]',[animal_data(:).SEM]',[animal_data(:).n],...
        [session_data(:).data]',[session_data(:).SEM]',[session_data(:).n]},NaN),...
        'VariableNames',cat(2,sdfvarnames{:}));
    writetable(T,sdfnames{f},'Delimiter','\t','WriteRowNames',false)    
end

% Print stats to console
fprintf('DREADD effects per animal:\n');
for nmds = 1:length(labels)
    fprintf('%s data vs. nov: p = %.4g\n',labels{nmds},pA(nmds));
end
fprintf('DREADD effects per session:\n');
for nmds = 1:length(labels)
    fprintf('%s data vs. nov: p = %.4g\n',labels{nmds},pS(nmds));
end

% SAVE FIGURE
pos = get(f1,'Position');
set(f1,'PaperPositionMode','Auto','PaperUnits','points','PaperSize',[pos(3), pos(4)])
print(f1,sprintf(['SOM-DREADD' '.pdf']),'-dpdf','-painters')

clear mds nmds names data cellinfo animal_data session_data...
    animalRes sessionRes parameter axistitle labels pA pS rdata vrange...
    ACTthr combination f fields isoutlier lims plotcats pPFthr pSIthr restable...
    selcats ylabels

%% K) Spatial Correlation PV:DREADD bsl vs clz
% SETTING PARAMETERS
ACTthr = 1/60;
pSIthr = .05;
lims = [-.5 1]; % Plotting range for each parameter. Needs to match fields

f1 = figure('position',[100 100 180 180]); hold on;

% FN Correlations baseline vs. cloz
[res,~,cellinfo] = metacorr_pairwise2(DREADDpv,'categories',1:4,'active',ACTthr,...
    'SI',pSIthr,'combination','or','plotresults',false);

data(1,:) = squeeze(res(1,2,:));
data(2,:) = squeeze(res(3,4,:));

% No outlier removal needed here - the range of correlation values is
% constrained enough!

% Plot by-animal and by-session data, get results in structures
[session_data, animal_data] = scatter_by_animal(...
    data, cellinfo, lut);

ylabel('clozapine');
xlabel('baseline');
%title('PV-Cre DREADD FN-correlation');

xlim(lims); ylim(lims);
plot(lims,lims,'-.','color',[.5 .5 .5]);
vrange = lims(2)-lims(1);

% Save source data
sdfvarnames{1} = strcat('PV-h4MDi animal',{' bsl',' clz',' bsl SEM',' clz SEM',' N cells'});
sdfvarnames{2} = strcat('PV-h4MDi session',{' bsl',' clz',' bsl SEM',' clz SEM',' N cells'});
T = array2table(catuneven({[animal_data(:).data]',[animal_data(:).SEM]',[animal_data(:).n],...
    [session_data(:).data]',[session_data(:).SEM]',[session_data(:).n]},NaN),...
    'VariableNames',cat(2,sdfvarnames{:}));
writetable(T,'SupplFig2k','Delimiter','\t','WriteRowNames',false)

% Statistics
fprintf('PV-Cre DREADD FN-correlation (animals):\n');
data = [animal_data(:).data];
pA = pairedSampleTest(data(1,:),data(2,:),'labels',{'baseline','clozapine'},'verbose',true);
t = text(lims(1)+0.05*vrange, lims(2)-0.05*vrange,...
    sprintf('p=%.4f, N=%i',pA,length(find([animal_data(:).n]>0))));
t.FontSize = 8;

printf('PV-Cre DREADD FN-correlation (sessions):\n');
data = [session_data(:).data];
pS = pairedSampleTest(data(1,:),data(2,:),'labels',{'baseline','clozapine'},'verbose',true);
t = text(lims(1)+0.4*vrange, lims(1)+0.075*vrange,...
    sprintf('p=%.4f, N=%i',pS,length(find([session_data(:).n]>0))));
t.FontSize = 8;

% SAVE FIGURE
pos = get(f1,'Position');
set(f1,'PaperPositionMode','Auto','PaperUnits','points','PaperSize',[pos(3), pos(4)])
print(f1,sprintf(['PV-DREADD corr' '.pdf']),'-dpdf','-painters')

clear mds nmds names data cellinfo animal_data session_data...
    animalRes sessionRes parameter axistitle labels pA pS rdata vrange...
    ACTthr combination f fields isoutlier lims plotcats pPFthr pSIthr restable...
    selcats ylabels corr res f1

%% P) Spatial Correlation SOM:DREADD bsl vs clz
% SETTING PARAMETERS
ACTthr = 1/60;
pSIthr = .05;
lims = [-.5 1]; % Plotting range for each parameter. Needs to match fields

f1 = figure('position',[100 100 180 180]); hold on;

% FN Correlations baseline vs. cloz
[res,~,cellinfo] = metacorr_pairwise2(DREADDsom,'categories',1:4,'active',ACTthr,...
    'SI',pSIthr,'combination','or','plotresults',false);

data(1,:) = squeeze(res(1,2,:));
data(2,:) = squeeze(res(3,4,:));

% No outlier removal needed here - the range of correlation values is
% constrained enough!

% Plot by-animal and by-session data, get results in structures
[session_data, animal_data] = scatter_by_animal(...
    data, cellinfo, lut);

ylabel('clozapine');
xlabel('baseline');
%title('SOM-Cre DREADD FN-correlation');

xlim(lims); ylim(lims);
plot(lims,lims,'-.','color',[.5 .5 .5]);
vrange = lims(2)-lims(1);

% Save source data
sdfvarnames{1} = strcat('SOM-h4MDi animal',{' bsl',' clz',' bsl SEM',' clz SEM',' N cells'});
sdfvarnames{2} = strcat('SOM-h4MDi session',{' bsl',' clz',' bsl SEM',' clz SEM',' N cells'});
T = array2table(catuneven({[animal_data(:).data]',[animal_data(:).SEM]',[animal_data(:).n],...
    [session_data(:).data]',[session_data(:).SEM]',[session_data(:).n]},NaN),...
    'VariableNames',cat(2,sdfvarnames{:}));
writetable(T,'SupplFig2p','Delimiter','\t','WriteRowNames',false)

% Statistics
fprintf('SOM-Cre DREADD FN-correlation (animals):\n');
data = [animal_data(:).data];
pA = pairedSampleTest(data(1,:),data(2,:),'labels',{'baseline','clozapine'},'verbose',true);
t = text(lims(1)+0.05*vrange, lims(2)-0.05*vrange,...
    sprintf('p=%.4f, N=%i',pA,size(data,2)));
t.FontSize = 8;

printf('SOM-Cre DREADD FN-correlation (sessions):\n');
data = [session_data(:).data];
pS = pairedSampleTest(data(1,:),data(2,:),'labels',{'baseline','clozapine'},'verbose',true);
t = text(lims(1)+0.4*vrange, lims(1)+0.075*vrange,...
    sprintf('p=%.4f, N=%i',pS,size(data,2)));
t.FontSize = 8;

% SAVE FIGURE
pos = get(f1,'Position');
set(f1,'PaperPositionMode','Auto','PaperUnits','points','PaperSize',[pos(3), pos(4)])
print(f1,sprintf(['SOM-DREADD corr' '.pdf']),'-dpdf','-painters')

clear mds nmds names data cellinfo animal_data session_data...
    animalRes sessionRes parameter axistitle labels pA pS rdata vrange...
    ACTthr combination f fields isoutlier lims plotcats pPFthr pSIthr restable...
    selcats ylabels corr res f1


%% SUPPLEMENTARY FIGURE 3
%% A,B,C) Running speed profiles SOM-Cre vs. PV-cre mice and # laps run
SOMds = {'CA1SOM','CA3SOM','DGSOM'};
PVds = {'CA1PV','CA3PV','DGPV'};
colors = {[1 0 1], [1 .5 0]};
clear PvVhist SomVhist PvVmean SomVmean
it = 0;

for nmds = 1:length(SOMds)
    eval(sprintf('mds = %s;',SOMds{nmds}))
    for nds = 1:length(mds)
        it = it+1;
        fr = mds{nds}.metadata.categories{1}.acquisition_rate(1);
        if iscell(fr); fr = fr{:}; end
        
        y = cat(2,mds{nds}.metadata.categories{1}.y{:},...
            mds{nds}.metadata.categories{2}.y{:})*200;
        v = movmedian(abs(diff(movmedian(y,round(fr/3))))*fr,fr);
        
        % Find lap ends
        y = movmedian(y,round(fr/3));
        ye = conv(y>prctile(y,90),true(1,round(fr)),'same');
        ys = conv(y<prctile(y,10),true(1,round(fr)),'same');
        
        %figure; hold on; plot(y); plot([0 diff(ye&ys)==1]*100,'r');
        
        SomVhist(it,:) = histcounts(v,1:31)/length(v);
        SomVmean(it) = nanmean(v);
        SomLaps(it) = length(find([0 diff(ye&ys)==1]));
    end
end

it = 0;
for nmds = 1:length(PVds)
    eval(sprintf('mds = %s;',PVds{nmds}))
    for nds = 1:length(mds)
        it = it+1;
        fr = mds{nds}.metadata.categories{1}.acquisition_rate(1);
        if iscell(fr); fr = fr{:}; end
        
        y = cat(2,mds{nds}.metadata.categories{1}.y{:},...
            mds{nds}.metadata.categories{2}.y{:})*200;
        v = movmedian(abs(diff(movmedian(y,round(fr/3))))*fr,fr);
        
        % Find lap starts and ends
        y = movmedian(y,round(fr/3));
        ye = conv(y>prctile(y,90),true(1,round(fr)),'same');
        ys = conv(y<prctile(y,10),true(1,round(fr)),'same');
        
        %figure; hold on; plot(y); plot([0 diff(ye&ys)==1]*100,'r');
        
        PvVhist(it,:) = histcounts(v,1:31)/length(v);
        PvVmean(it) = nanmean(v);
        PvLaps(it) = length(find([0 diff(ye&ys)==1]));
    end
end

figure; hold on;
errorbar(mean(SomVhist,1)*100,...
    std(SomVhist,[],1)./sqrt(size(SomVhist,1))*100,'o-','Color',colors{1});
errorbar(mean(PvVhist,1)*100,...
    std(PvVhist,[],1)./sqrt(size(PvVhist,1))*100,'o-','Color',colors{2});
xlabel('Speed (cm/s)');
ylabel('Fraction of time (%)');
ylim([0 4]);
legend({'SOM','PV'});


% Save source data file
T = table([mean(SomVhist,1)*100]',[std(SomVhist,[],1)./sqrt(size(SomVhist,1))*100]',...
    [mean(PvVhist,1)*100]',[std(PvVhist,[],1)./sqrt(size(PvVhist,1))*100]',...
    'VariableNames',{'SOM Speed (cm/s)','SOM speed SEM','PV Speed (cm/s)','PV speed SEM'});
writetable(T,'SupplFig3a','Delimiter','\t');

% Mean speed comparison
figure('position',[200 200 150 220]); hold on;
scatter(randn(length(SomVmean),1)/10+1,SomVmean,[],colors{1},'.')
scatter(randn(length(PvVmean),1)/10+2,PvVmean,[],colors{2},'.')
errorbar(1, mean(SomVmean,2), std(SomVmean,[],2)./sqrt(size(SomVmean,2)),'o-','Color','k');
errorbar(2, mean(PvVmean,2), std(PvVmean,[],2)./sqrt(size(PvVmean,2)),'o-','Color','k');
xlim([0.5 2.5]);
xticks([1 2]);
xticklabels({'SOM','PV'});
ylabel('Mean running speed (cm/s)');

% Save source data file
T = array2table(catuneven({SomVmean,PvVmean},NaN),'VariableNames',{'SOM','PV'});
writetable(T,'SupplFig3b','Delimiter','\t');

% Number of laps completed
figure('position',[400 200 150 220]); hold on;
scatter(randn(length(SomLaps),1)/10+1,SomLaps,[],colors{1},'.')
scatter(randn(length(PvLaps),1)/10+2,PvLaps,[],colors{2},'.')
errorbar(1, mean(SomLaps,2), std(SomLaps,[],2)./sqrt(size(SomLaps,2)),'o-','Color','k');
errorbar(2, mean(PvLaps,2), std(PvLaps,[],2)./sqrt(size(PvLaps,2)),'o-','Color','k');
xlim([0.5 2.5]);
xticks([1 2]);
xticklabels({'SOM','PV'});
ylabel('Total # laps / session');

% Save source data file
T = array2table(catuneven({SomLaps,PvLaps},NaN),'VariableNames',{'SOM','PV'});
writetable(T,'SupplFig3c','Delimiter','\t');

fprintf('Comparing mean running speeds:\n');
unpairedSampleTest(SomVmean,PvVmean,'verbose',true,'labels',{'SOM','PV'});

fprintf('Comparing number of laps per experiment:\n');
unpairedSampleTest(SomLaps,PvLaps,'verbose',true,'labels',{'SOM','PV'});

clear v y names nmds mds nds it PvVhist SomVhist PvVmean SomVmean

%% D) Running speed profile Fam-Nov transition
names = {'CA1PV','CA1SOM','CA3PV','CA3SOM','DGPV','DGSOM'};
intv = 30; % Around transition, in seconds
fts = -intv:.1:-.1; % Timestamps, standartized
nts = 0:.1:intv; % Timestamps, standartized

allv = []; it = 0;

for nmds = 1:length(names)
    eval(sprintf('mds = %s;',names{nmds}))
    for nds = 1:length(mds)
        for r = 1:length(mds{nds}.metadata.categories{1}.y)
            fr = mds{nds}.metadata.categories{1}.acquisition_rate(r);
            
            % Get speed for fam and nov
            vf = mds{nds}.metadata.categories{1}.y{r}*200;
            vf = abs(diff(vf))*fr;
            %vf = movmedian(abs(diff(movmedian(vf,round(fr/3))))*fr,.25);
            vn = mds{nds}.metadata.categories{2}.y{r}*200;
            vn = abs(diff(vn))*fr;
            %vn = movmedian(abs(diff(movmedian(vn,round(fr/3))))*fr,.25);
            
            if length(vf)<intv*fr+fr | length(vn)<=intv*fr+fr
                continue
            end
            
            % Clip speed traces to interval
            vf = vf(end-round(intv*fr):end);
            vf = interp1(linspace(-intv,-1/fr,length(vf)),vf,fts);
            vn = vn(1:round(intv*fr));
            vn = interp1(linspace(1/fr,intv,length(vn)),vn,nts);
            
            % Remove artificial 'spikes' in speed related to enviroment
            % switching

            it = it+1;
            allv(it,:) = cat(2,vf,vn);
            
            % Remove artificial 'spikes' in speed related to pos reset
            allv(it,allv(it,:)>100) = NaN;
        end
    end
end

figure('Position',[400 400 400 400]); hold on;
shadedErrorBar(cat(2,fts,nts),nanmean(allv,1),nanstd(allv,[],1)./sqrt(size(allv,1)));
xlabel('Time from switch (s)');
ylabel('Speed (cm/s)');

figure;
imagesc(cat(2,fts,nts),1:size(allv,2),allv);
caxis([0 50]);
xlabel('Time from switch (s)');
ylabel('Run #');

clear v c it fr fts intv mds names nds nmds nts r vf vn

%% E) Running speed over trial number
names = {'CA1PV','CA1SOM','CA3PV','CA3SOM','DGPV','DGSOM'};
labels = {'F1','N1','F2','N2','F3','N3','F4','N4','F5','N5'};
allv = []; it = 0;

for nmds = 1:length(names)
    eval(sprintf('mds = %s;',names{nmds}))
    for nds = 1:length(mds)
        it = it+1;
        ep = 0;
        for r = 1:length(mds{nds}.metadata.categories{1}.y)
            for c = 1:2
                ep = ep+1;
                fr = mds{nds}.metadata.categories{c}.acquisition_rate(r);
                v = cat(2,mds{nds}.metadata.categories{c}.y{r})*200;
                v = movmedian(abs(diff(movmedian(v,round(fr/3))))*fr,fr);
                allv(it,ep) = nanmean(v);
                %vhist{nmds,nds,c} = histcounts(v,1:31)/length(v);
            end
        end
    end
end

% Plotting
figure('Position',[400 400 400 400]); hold on;
%plot(allv(:,1:10)','Color',[.5 .5 .5]);
errorbar(1:10, mean(allv(:,1:10),1), std(allv(:,1:10),1)./sqrt(size(allv,1)),...
    'o-','Color','k','CapSize',0);
xlim([.5 10.5]);
xticks(1:10);
xticklabels(labels);
xlabel('Run Epoch');
ylabel('Mean running speed (cm/s)');
ylim([0 20]);

% Save Source Data file
T = table([mean(allv(:,1:10),1)]',[std(allv(:,1:10),1)./sqrt(size(allv,1))]',...
    'VariableNames',{'Speed (cm/s)','Speed SEM'},'RowNames',labels);
writetable(T,'SupplFig3e','Delimiter','\t','WriteRowNames',true);

% Statistics
kruskalwallisNdunns(allv(:,1:10),labels)

% Pearson correlation
[R,p] = corrcoef(repmat(1:10,54,1),allv(:,1:10))

clear v c it ep fr mds names nds nmds r labels

%% F) Running Speed profiles across all DS, mean +/- SEM
names = {'CA1PV','CA1SOM','CA3PV','CA3SOM','DGPV','DGSOM'};
it = 0;
figure; hold on;
for nmds = 1:length(names)
    eval(sprintf('mds = %s;',names{nmds}))
    for nds = 1:length(mds)
        it = it+1;
        for c = 1:2
            fr = mds{nds}.metadata.categories{c}.acquisition_rate(1);
            v = cat(2,mds{nds}.metadata.categories{c}.y{:})*200;
            v = movmedian(abs(diff(movmedian(v,round(fr/3))))*fr,fr);
            Vhist{it,c} = histcounts(v,1:31)/length(v);
        end
    end
end
errorbar(mean(cat(1,Vhist{:,1}))*100,...
    std(cat(1,Vhist{:,1}))./sqrt(size(Vhist,1))*100,'o-','Color',[1 0 0]);
errorbar(mean(cat(1,Vhist{:,2}))*100,...
    std(cat(1,Vhist{:,2}))./sqrt(size(Vhist,1))*100,'o-','Color',[0 0 1]);
xlabel('Speed (cm/s)');
ylabel('Fraction of time (%)');
legend({'fam','nov'});

% Save Source Data file
T = table([mean(cat(1,Vhist{:,1}))*100]',[std(cat(1,Vhist{:,1}))./sqrt(size(Vhist,1))*100]',...
    [mean(cat(1,Vhist{:,2}))*100]',[std(cat(1,Vhist{:,2}))./sqrt(size(Vhist,1))*100]',...
    'VariableNames',{'Fam Speed (cm/s)','Fam speed SEM','Nov Speed (cm/s)','Nov speed SEM'});
writetable(T,'SupplFig3f','Delimiter','\t');

clear v names nmds mds nds it

%% G) Mean Running Speed fam vs. nov across all datasets
names = {'CA1PV','CA1SOM','CA3PV','CA3SOM','DGPV','DGSOM'};
allv = []; it = 0;

for nmds = 1:length(names)
    eval(sprintf('mds = %s;',names{nmds}))
    for nds = 1:length(mds)
        it = it+1;
        for c = 1:2
            fr = mds{nds}.metadata.categories{c}.acquisition_rate(1);
            v = cat(2,mds{nds}.metadata.categories{c}.y{:})*200;
            v = movmedian(abs(diff(movmedian(v,round(fr/3))))*fr,fr);
            allv(it,c) = nanmean(v);
            %vhist{nmds,nds,c} = histcounts(v,1:31)/length(v);
        end
    end
end

figure('Position',[400 400 110 250]); hold on;
plot(allv','Color',[.5 .5 .5]);
errorbar(1:2, mean(allv,1), std(allv,1)./sqrt(size(allv,1)),'o-','Color','k');
xlim([.75 2.25]);
xticks([1 2]);
xticklabels({'familiar','novel'});
ylabel('Mean running speed (cm/s)');
title(sprintf('p = %.2d',pairedSampleTest(allv(:,1),allv(:,2))));

% Statistics
pairedSampleTest(allv(:,1),allv(:,2),'verbose',true);

% Save Source Data file
T = table(allv(:,1),allv(:,2),'VariableNames',{'Fam mean Speed (cm/s)','Nov mean Speed (cm/s)'});
writetable(T,'SupplFig3g','Delimiter','\t');

clear v c it fr mds names nds mnmds

%% H) Activity differences fam vs. nov violin plots w/o shuffle
% See Fig 4A,B (code above).


%% SUPPLEMENTARY FIGURE 4
%% A) Example cells with circular tuning vectors
nsamp = length(CA1{1}.cells{1}.categories{1}.dFoY{1}); %Assuming identical spatial bins
xrad = ((1:nsamp)*2*pi/nsamp)';
xv = deg2rad(1:360);
colmap = jet;
%colmap(1,:) = [1 1 1];

names = {'CA1','CA1PV','CA1SOM','CA3PV','CA3SOM','DGPV','DGSOM'};
labels = {'CA1-PYR','CA1-PV','CA1-SOM','CA2/3-PV','CA2/3-SOM','DG-PV','DG-SOM'};
smp = [4,25; 5,2; 9,3; 4,13; 1,44; 6,3; 4,1];
colors = {[1 1 0],[0 .5 0],[.5 1 .5],[0 0 1],[.5 .5 1],[1 0 0],[1 .5 .5]};

f = figure;
set(f,'Units','Inch','Position',[2 2 8.1 2.6],'defaultAxesFontSize',6);

for nmds = 1:length(names)
    eval(sprintf('mds = %s;',names{nmds}))
   
    thisDFoY = nanmean(cat(2,mds{smp(nmds,1)}.cells{smp(nmds,2)}.categories{1}.dFoY{:}),2);
    thisA = mds{smp(nmds,1)}.cells{smp(nmds,2)}.tva;
    thisTVL = mds{smp(nmds,1)}.cells{smp(nmds,2)}.tvl;
    
    % Linear maps across runs
    subplot(2,7,nmds);
    dFoY = cat(2,mds{smp(nmds,1)}.cells{smp(nmds,2)}.categories{1}.dFoY{:})';
    dFoY = movmean(dFoY,2,2);
    xl = linspace(0,4,size(dFoY,2));
    yl = 1:size(dFoY,1);
    imagesc(xl,yl,dFoY);
    xlabel('Track distance (m)');
    %ylabel('Run #');
    colormap(colmap);
    caxis([-.2 .8]);    
    title(labels{nmds})
    
    % Plot radial distribution and vector
    subplot(2,7,nmds+7);
    polarplot(xrad, thisDFoY,'Color',colors{nmds},'LineWidth',1.4);
    hold on;

    yv = zeros(1,360);
    yv(round(mod(rad2deg(thisA),360))) = thisTVL;
    polarplot(xv,yv,'Color','k','LineWidth',2);
    
    Ax = gca; % current axes
    Ax.ThetaTickLabel = [];
end

pos = get(f,'Position');
set(f,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(f,'Circular_tuning_examples.pdf','-dpdf','-r0')

clear mds nmds xv yv xl yl thisA thisTVL thisDFoY colors xrad names smp nsamp labesl pos f dFoY Ax colmap

%% B) Circular tuning vs. shuffle
names = {'CA1PV','CA1SOM','CA3PV','CA3SOM','DGPV','DGSOM'};
labels = ({'CA1 PV','CA1 SOM','CA2/3 PV','CA2/3 SOM','DG PV','DG SOM'});
param = 'tvl';%'MovImmRatio','speedmod','trialvar','sessionstab','spatial_coherence','spatialinfo','SIperAUC','tvl','positional_info'
rparam = sprintf([param '_r']);
ylims = [-.02 .62];
yl = 'Positional Information';%'Tuning-vector length';
OutlierThr = 3; % In interquartile ranges, use tukey procedure
colors = {[0 .5 0],[.5 1 .5],[0 0 1],[.5 .5 1],[1 0 0],[1 .5 .5]};
randcolor = [.8 .8 .5];
SI = 1;

clear values rvalues
for nmds = 1:length(names)
    eval(sprintf('mds = %s;',names{nmds}))
    clear val
    
    val(:,1) = metastatistics(mds,'cell',param,'categories',1,...
        'rescategories',1,'plotting',false,'SI',SI)';
    val(:,2) = metastatistics(mds,'cell',rparam,'categories',1,...
        'rescategories',1,'plotting',false)';
    
    val = tukeyOutlierRemoval(val,OutlierThr);
    
    p(nmds) = signrank(val(:,1),val(:,2));
    N(nmds) = size(val,1);
    
    values{nmds} = val(:,1);
    rvalues{nmds} = val(:,2);
end
close all

% Build 'global' array with all groups for stakeplot -- note, here the
% 'real' values are first, bootstrap second
valA = cat(1,values,rvalues);
valA = catuneven(reshape(valA,1,[]),NaN);

% Save source data files
sdfnames = cat(1,names,strcat(names,' rand'));
sdfnames = reshape(sdfnames,1,[]);
T = array2table(valA,'VariableNames',sdfnames);
writetable(T,'SupplFig4b','Delimiter','\t','WriteRowNames',false)  

f3 = figure('Units','centimeters','position',[4,4,6,4.5],'defaultAxesFontSize',8);
%title(param);
stakeplot(valA,'color',[.3 .3 .3],'LineAlpha',.2,'spread',.45,'jitter',0,'LineWidth',.5);
distributionPlot(values,'histOri','left','color',colors,'widthDiv',[2 1],'showMM',6,'FaceAlpha',.7);
distributionPlot(rvalues,'histOri','right','color',randcolor,'widthDiv',[2 2],'showMM',6,'FaceAlpha',.7);
%xticks([])
xticklabels(labels);
xtickangle(45);
ylim(ylims);
ylabel(yl);
%set(gca, 'YScale', 'log')

% Create PDF output
% set(findobj(f3,'type','text'),'fontsize',7)
pos = get(f3,'Position');
set(f3,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(f3,sprintf('%s_vs_chance.pdf',param),'-dpdf','-painters')

% Plot KruskalWallisAnova across grops to console
kruskalwallisNdunns(catuneven(values,NaN),names);

% Print group-wise N and p-values to console
fprintf('Number of observations (N)\n');
for nmds = 1:length(names)
    fprintf('%s: %i\n',names{nmds},N(nmds));
end

fprintf('\nSigned ranks p-value\n');
for nmds = 1:length(names)
    fprintf('%s: %d\n',names{nmds},p(nmds));
end

clear param rparam rval rvalues nmds mds yl names val...
    colors oval f3 pos SI ylims N p OutlierThr randcolor values

%% C) Positional info vs. shuffle
names = {'CA1PV','CA1SOM','CA3PV','CA3SOM','DGPV','DGSOM'};
labels = ({'CA1 PV','CA1 SOM','CA2/3 PV','CA2/3 SOM','DG PV','DG SOM'});
param = 'positional_info';%'MovImmRatio','speedmod','trialvar','sessionstab','spatial_coherence','spatialinfo','SIperAUC','tvl','positional_info'
rparam = sprintf([param '_r']);
ylims = [-.02 .62];
yl = 'Positional Information';%'Tuning-vector length';
OutlierThr = 3; % In interquartile ranges, use tukey procedure
colors = {[0 .5 0],[.5 1 .5],[0 0 1],[.5 .5 1],[1 0 0],[1 .5 .5]};
randcolor = [.8 .8 .5];
SI = 1;

clear values rvalues
for nmds = 1:length(names)
    eval(sprintf('mds = %s;',names{nmds}))
    clear val
    
    val(:,1) = metastatistics(mds,'cell',param,'categories',1,...
        'rescategories',1,'plotting',false,'SI',SI)';
    val(:,2) = metastatistics(mds,'cell',rparam,'categories',1,...
        'rescategories',1,'plotting',false)';
    
    val = tukeyOutlierRemoval(val,OutlierThr);
    
    p(nmds) = signrank(val(:,1),val(:,2));
    N(nmds) = size(val,1);
    
    values{nmds} = val(:,1);
    rvalues{nmds} = val(:,2);
end
close all

% Build 'global' array with all groups for stakeplot
valA = cat(1,values,rvalues);
valA = catuneven(reshape(valA,1,[]),NaN);

% Save source data files
sdfnames = cat(1,names,strcat(names,' rand'));
sdfnames = reshape(sdfnames,1,[]);
T = array2table(valA,'VariableNames',sdfnames);
writetable(T,'SupplFig4c','Delimiter','\t','WriteRowNames',false)  

f3 = figure('Units','centimeters','position',[4,4,6,4.5],'defaultAxesFontSize',8);
%title(param);
stakeplot(valA,'color',[.3 .3 .3],'LineAlpha',.2,'spread',.45,'jitter',0,'LineWidth',.5);
distributionPlot(values,'histOri','left','color',colors,'widthDiv',[2 1],'showMM',6,'FaceAlpha',.7);
distributionPlot(rvalues,'histOri','right','color',randcolor,'widthDiv',[2 2],'showMM',6,'FaceAlpha',.7);
%xticks([])
xticklabels(labels);
xtickangle(45);
ylim(ylims);
ylabel(yl);
%set(gca, 'YScale', 'log')

% Create PDF output
% set(findobj(f3,'type','text'),'fontsize',7)
pos = get(f3,'Position');
set(f3,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(f3,sprintf('%s_vs_chance.pdf',param),'-dpdf','-painters')

% Plot KruskalWallisAnova across grops to console
kruskalwallisNdunns(catuneven(values,NaN),names);

% Print group-wise N and p-values to console
fprintf('Number of observations (N)\n');
for nmds = 1:length(names)
    fprintf('%s: %i\n',names{nmds},N(nmds));
end

fprintf('\nSigned ranks p-value\n');
for nmds = 1:length(names)
    fprintf('%s: %d\n',names{nmds},p(nmds));
end

clear param rparam rval rvalues nmds mds yl names val...
    colors oval f3 pos SI ylims N p OutlierThr randcolor values

%% D-G) Spatial tuning properties vs. Principal cells
selected = 6; % Select parameters and sdfnames
param = {'SIperAUC','tvl','positional_info','trialvar','sessionstab','sessionstab'};
sdfname = {'SupplFig4d','SupplFig4e','SupplFig4f','SupplFig4g','SupplFig5c','SupplFig5f'};
yl = {'Spatial information (normalized)','Tuning-vector length','Positional information',...
    'Trial-to-trial correlation','Within session stability','Within session stability'};
SI = 1; if selected == 6; SI = 0.05; end
names = {'CA1','CA1PV','CA1SOM','CA3PV','CA3SOM','DGPV','DGSOM'};
labels = {'CA1-PYR','CA1-PV','CA1-SOM','CA2/3-PV','CA2/3-SOM','DG-PV','DG-SOM'};
%param = 'sessionstab';%'MovImmRatio','speedmod','trialvar','sessionstab','spatial_coherence','tvl'
ylims = {[-.2 4; -.2 2],[0 2; 0 .3],[0 .6; 0 .6],[-.5 1; -.25 .5],[-.5 1; -.5 1],[-.5 1; -.5 1]};
%ylims = [-.02 .6; -.02 .6];

param = param{selected};
yl = yl{selected};
ylims = ylims{selected};

OutlierThr = 3; % In interquartile ranges, use tukey procedure
colors = {[1 1 0],[0 .5 0],[.5 1 .5],[0 0 1],[.5 .5 1],[1 0 0],[1 .5 .5]};

for nmds = 1:length(names)
    eval(sprintf('mds = %s;',names{nmds}))
    val = metastatistics(mds,'cell',param,'categories',1,...
        'rescategories',1,'plotting',false,'SI',SI)';
    
    % outlier removal
    val = tukeyOutlierRemoval(val',OutlierThr);
    
    values{nmds} = val;
    N(nmds) = length(find(~isnan(val)));
end

f3 = figure;
set(f3,'Units','centimeters','position',[5 5 6.7 5.0],'defaultAxesFontSize',6);
beeswarmplot(catuneven(values,NaN),'color',[.7 .7 .7],'markersize',4,'MarkerFaceAlpha',.5);
distributionPlot(values,'color',colors,'showMM',6,'FaceAlpha',.5);
xticklabels(labels);
xtickangle(45);
ylabel(yl);
ylim(ylims(1,:));
set(gca,'FontSize',6);
%set(gca, 'YScale', 'log')

pos = get(f3,'Position');
set(f3,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(f3,sprintf('%s-SI%.2f.pdf',param,SI),'-dpdf','-r0')

% Save source data file
T = array2table(catuneven(values,NaN),'VariableNames',names);
writetable(T,sdfname{selected},'Delimiter','\t','WriteRowNames',false)

% Print inset of IN values only
f4 = figure;
set(f4,'defaultAxesFontSize',6,'Units','Inches','Position',[2 2 1.5 1.3]);
beeswarmplot(catuneven(values(2:7),NaN),'color',[.7 .7 .7],'markersize',4,'MarkerFaceAlpha',.5);
distributionPlot(values(2:7),'color',colors(2:7),'showMM',6,'FaceAlpha',.5);
xticklabels(labels(2:7));
xtickangle(45);
ylim(ylims(2,:));
set(gca,'FontSize',6);
%set(gca, 'YScale', 'log')

pos = get(f4,'Position');
set(f4,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(f4,sprintf('%s-SI%.2f-inset.pdf',param,SI),'-dpdf','-r0')

% Print group-wise N and p-values to console
fprintf('Number of observations (N)\n');
for nmds = 1:length(names)
    fprintf('%s: %i\n',names{nmds},N(nmds));
end

kruskalwallisNdunns(catuneven(values,NaN),names);

clear param rparam rval rvalues nmds mds yl names val Foutliers...
    colors oval values OutlierThr pos colors sdfname

%% REVIEWER ONLY: "Negative spatial tuning"
% Get 'place fields' as points significantly deviating from the
% distribution of fluorescence across all positions (at least 3 consecutive) and
% assess the fraction of cells having 'positive', 'negative' or 'both'
% fields.
names = {'CA1PV','CA1SOM','CA3PV','CA3SOM','DGPV','DGSOM'};
labels = {'CA1-PV','CA1-SOM','CA2/3-PV','CA2/3-SOM','DG-PV','DG-SOM'};
f1 = figure;
set(f1,'Units','Inches','Position',[2 1 6 1.2],'defaultAxesFontSize',12);

for nmds = 1:length(names)
    eval(sprintf('mds = %s;',names{nmds}))
    nPosFields = metastatistics(mds,'cell','nPosFields','categories',1,...
        'rescategories',1,'plotting',false)';
    nNegFields = metastatistics(mds,'cell','nNegFields','categories',1,...
        'rescategories',1,'plotting',false)';
    hasBoth = nPosFields>0 & nNegFields>0;
    hasPos = nPosFields>0 & nNegFields==0;
    hasNeg = nPosFields==0 & nNegFields>0;
    %hasNone = nPosFields == 0 & nNegFields == 0;
    
    tbl(1:3,nmds) = cat(1,length(find(hasPos)), length(find(hasNeg)),length(find(hasBoth)));
    plbs = tbl(1:3,nmds)/sum(tbl(1:3,nmds))*100;
    plbs = {sprintf('%.1f%%',plbs(1)), sprintf('%.1f%%',plbs(2)),...
        sprintf('%.1f%%',plbs(3))};
    
    figure(f1);
    subplot(1,length(names),nmds)
    pie(tbl(:,nmds),plbs);
    ax = gca;
    ax.Colormap = [1 0 0; 0 0 1; .5 .5 .5];
    title(labels{nmds},'FontSize',8);
    clear plbs
end

set(findobj(f1,'type','text'),'fontsize',7)

% Save as pdf
pos = get(f1,'Position');
set(f1,'PaperPositionMode','Auto','PaperUnits','Inches',...
    'PaperSize',[pos(3), pos(4)],'Renderer','Painters')
print(f1,'PosNegPFfractions.pdf','-dpdf','-r0')

[Chi2,p] = chi2tbltest(tbl);
fprintf('Chi2 = %.2f; p = %d\n',Chi2,p);
clear names speed_P speed_R tuned_Pos tuned_Neg untuned tbl nmds mds ax Chi2 p f1 labels


%% SUPPLEMENTARY FIGURE 5 
%% A,B,D,E) Placemaps first-second half split + correlation matrix
%names = {'CA1PV','CA1SOM','CA3PV','CA3SOM','DGPV','DGSOM'};
SI = 1; % Set to 0.05 for D,E
names = {'CA1'};
f1 = figure('Position', [100 100 700 500]);
fc = figure('Position', [100 100 330 500]);

for nmds = 1:length(names)
    eval(sprintf('mds = %s;',names{nmds}))
    
    for nds = 1:length(mds)
        % Duplicate dataset
        mdss{nds} = mds{nds};
        ntraces = length(mds{nds}.cells{1}.categories{1}.dFoY);
        mdss{nds}.metadata.categories{2} = mds{nds}.metadata.categories{1};
        mdss{nds}.metadata.categories{1}.x = ...
            mds{nds}.metadata.categories{1}.x(1:round(ntraces/2));
        mdss{nds}.metadata.categories{2}.x = ...
            mds{nds}.metadata.categories{1}.x(round(ntraces/2)+1:end);
        mdss{nds}.metadata.categories{1}.y = ...
            mds{nds}.metadata.categories{1}.y(1:round(ntraces/2));
        mdss{nds}.metadata.categories{2}.y = ...
            mds{nds}.metadata.categories{1}.y(round(ntraces/2)+1:end);
        mdss{nds}.metadata.categories{1}.ft = ...
            mds{nds}.metadata.categories{1}.ft(1:round(ntraces/2));
        mdss{nds}.metadata.categories{2}.ft = ...
            mds{nds}.metadata.categories{1}.ft(round(ntraces/2)+1:end);
        
        for n = 1:length(mds{nds}.cells)
            mdss{nds}.cells{n}.transientrate(2) = 0;
            mdss{nds}.cells{n}.Placefield_P(2) = 1;
            mdss{nds}.cells{n}.spatial_P(2) = mdss{nds}.cells{n}.spatial_P(1);
            
            mdss{nds}.cells{n}.categories{1}.dFoT = ...
                mds{nds}.cells{n}.categories{1}.dFoT(1:round(ntraces/2));
            mdss{nds}.cells{n}.categories{2}.dFoT = ...
                mds{nds}.cells{n}.categories{1}.dFoT(round(ntraces/2)+1:end);
            mdss{nds}.cells{n}.categories{1}.dFoY = ...
                mds{nds}.cells{n}.categories{1}.dFoY(1:round(ntraces/2));
            mdss{nds}.cells{n}.categories{2}.dFoY = ...
                mds{nds}.cells{n}.categories{1}.dFoY(round(ntraces/2)+1:end);
            mdss{nds}.cells{n}.categories{1}.zscored = ...
                mds{nds}.cells{n}.categories{1}.zscored(1:round(ntraces/2));
            mdss{nds}.cells{n}.categories{2}.zscored = ...
                mds{nds}.cells{n}.categories{1}.zscored(round(ntraces/2)+1:end);
        end
    end
    [dFoYplot, sidcs] = plot_placemaps_meta(mdss,1,'categories',1,...
        'signal','zscored','bins',0.1:.025:1.9,'runmean',5,'SI',SI);
    %dFoYplot = movmean(dFoYplot,3,2);
    
    figure(f1)
    xv = (.1:.025:1.9)*2; % In meters
    yv = 1:length(sidcs);
    
    subplot(3,4,2*(nmds-1) + 1)
    imagesc(xv,yv,dFoYplot(:,sidcs,1)');
    xlabel('Track distance (m)'); ylabel('Cell #');
    caxis([-.1 2]); colormap(hot)
    title(sprintf('%s - 1st half',names{nmds}));
    
    subplot(3,4,2*(nmds-1) + 2)
    imagesc(xv,yv,dFoYplot(:,sidcs,2)');
    xlabel('Track distance (m)'); ylabel('Cell #');
    caxis([-.1 2]); colormap(hot)
    title(sprintf('%s - 2nd half',names{nmds}));
    %plot_placemaps_meta(GFPs,1,'categories',2,'signal','zscored','bins',0.175:0.025:2.0);
    
    figure(fc);
    % Z-score maps across space to account for rate-differences
    dFoYplot(:,:,1) = zscore(dFoYplot(:,:,1));
    dFoYplot(:,:,2) = zscore(dFoYplot(:,:,2));
    
    corrim = corr(dFoYplot(:,:,1)',dFoYplot(:,:,2)','rows','complete');
    subplot(3,2,nmds)
    imagesc(xv,xv,corrim);
    xlabel('Location 1st half (m)');
    ylabel('Location 2nd half (m)');
    caxis([-1 1])
    colormap(jet)
    title(names{nmds});
end

clear ds n ntraces xv yv nds nmds mds mdss names order sidcs dFoYplot f1

%% C,F) Within session stability
% See Supplementary 4, D-G) Set 'sessionstab' and 'SI=1' or 'SI=0.05' for
% C, and F, respectively.


%% SUPPLEMENTARY FIGURE 6
%% B) Individual example of ML SOM axon, multiple transitions
plot_FNtransition(AxonSomML{3},13,'dFoT_fneu',15);
xlim([-60 60]);

%% C) All ML SOM axons, FN transition w/o duplicates
smoothFct = 5; % For duplicate removal
corrThr = 0.7; % For duplicate removal
signaltype = 'dFoT_fneu'; % 'dFoT','dFoT_fneu'

% Get Hilus axon data
it = 0;
for ds = 1:length(AxonSomML)
    % Remove highly correlated boutons on local copy of the dataset
    MlWoD{ds} = remove_duplicates(AxonSomML{ds},corrThr,...
        'categories',[1 2],'signal',signaltype,'smoothFct',smoothFct);
end

[psth,fnratio,ts] = FNtransition_PSTH(MlWoD,signaltype);
[~,ord] = sort(fnratio);
figure; imagesc(ts,1:size(psth,2),psth(:,ord)'); 
caxis([0 1]); xlim([-60 60]); colormap('hot');
title('ML SOM axons');
ylabel('Axon #');
xlabel('Time (s)');

clear MlWoD ts psth ord fnratio smoothFct corrThr signaltype

%% E) Individual example of hilus SOM axon, multiple transitions
% Good expls: 1,18; 2,9 2,34; 5, 58 5, 81, 5,103, 
plot_FNtransition(AxonSomHilus{5},103,'dFoT_fneu',8);
xlim([-60 60]);

%% F) All hilar SOM axons, FN transition w/o duplicates
smoothFct = 5; % For duplicate removal
corrThr = 0.7; % For duplicate removal
signaltype = 'dFoT_fneu'; % 'dFoT','dFoT_fneu'

% Get Hilus axon data
for ds = 1:length(AxonSomHilus)
    % Remove highly correlated boutons on local copy of the dataset
    HilWoD{ds} = remove_duplicates(AxonSomHilus{ds},corrThr,...
        'categories',[1 2],'signal',signaltype,'smoothFct',smoothFct);
end

[psth,fnratio,ts] = FNtransition_PSTH(HilWoD,signaltype);
[~,ord] = sort(fnratio);
figure; imagesc(ts,1:size(psth,2),psth(:,ord)'); 
caxis([0 1]); xlim([-60 60]); colormap('hot');
title('Hilar SOM axons');
ylabel('Axon #');
xlabel('Time (s)');

clear ts ord fnratio smoothFct corrThr signaltype HilWoD psth

%% G) Head-on comparison of FN effect Hilus vs. ML axons
outlier_thr = 3; % 3 extreme, 1.5 moderate, Inf to disable
smoothFct = 5; % For duplicate removal
corrThr = 0.7; % For duplicate removal
ylims = [0 2.5];
colors = {[1 .7 .7],[1 .85 .85]};
labels = {'Hilus','ML'};
signaltype = 'dFoT_fneu'; % 'dFoT','dFoT_fneu'
axislabel = 'Activity Ratio (nov/fam)';

% Get Hilus axon data
it = 0;
for ds = 1:length(AxonSomHilus)
    % Remove highly correlated boutons on local copy of the dataset
    HilWoD{ds} = remove_duplicates(AxonSomHilus{ds},corrThr,...
        'categories',[1 2],'signal',signaltype,'smoothFct',smoothFct);
    
    % Iterate through unique axons, get mean activity in FAM and NOV.
    for n = 1:length(HilWoD{ds}.cells)
        it = it+1;
        for c = [1 2] % Fam and Nov only
            signal = cat(2,HilWoD{ds}.cells{n}.categories{c}.(signaltype){:});
            moving = cat(2,HilWoD{ds}.metadata.categories{c}.moving{:});
            HilusDF(it,c) = nanmean(signal(moving));
        end
    end
end

% Get ML axon data
it = 0;
for ds = 1:length(AxonSomML)
    % Remove highly correlated boutons on local copy of the dataset
    MlWoD{ds} = remove_duplicates(AxonSomML{ds},corrThr,...
        'categories',[1 2],'signal',signaltype,'smoothFct',smoothFct);
    
    % Iterate through unique axons, get mean activity in FAM and NOV.
    for n = 1:length(MlWoD{ds}.cells)
        it = it+1;
        for c = [1 2] % Fam and Nov only
            signal = cat(2,MlWoD{ds}.cells{n}.categories{c}.(signaltype){:});
            moving = cat(2,MlWoD{ds}.metadata.categories{c}.moving{:});
            MlDF(it,c) = nanmean(signal(moving));
        end
    end
end

HilusRatio = HilusDF(:,2)./HilusDF(:,1); % Increase in novel > 1
MlRatio = MlDF(:,2)./MlDF(:,1);

HilusRatio = tukeyOutlierRemoval(HilusRatio,outlier_thr);
MlRatio = tukeyOutlierRemoval(MlRatio,outlier_thr);

% Save source data
T = array2table(catuneven({HilusRatio,MlRatio},NaN),'VariableNames',{'Hilus Ratio','ML ratio'});
writetable(T,'SupplFig6g','Delimiter','\t','WriteRowNames',false)

% % PLOTTING
f1 = figure; hold on;
set(f1,'Units','centimeters','defaultAxesFontSize',10,'Position',[4,4,5.6,5.6]);

distributionPlot({HilusRatio,MlRatio},'color',colors,'showMM',6,'distWidth',.9,...
    'addSpread',0,'histOpt',1,'divFactor',3,'FaceAlpha',.5);

plotSpread({HilusRatio,MlRatio},'spreadFcn',{'lin',100},'spreadWidth',.9,...
    'distributionColors',{[max(colors{1}-.2,0) .01],[max(colors{2}-.2,0) .01]},...
    'distributionMarkerSize',1,'binWidth',.05);

ylim(ylims);
ylabel(axislabel,'FontSize',10);
xticklabels(labels);
xtickangle(45);
xlim([.5 2.5])
set(gca,'FontSize',10);

% SAVE FIGURE
pos = get(f1,'Position');
set(f1,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[pos(3), pos(4)])
print(f1,sprintf('SOM ML vs. Hilus.pdf'),'-dpdf','-painters','-bestfit')

fprintf('Activity ratios (nov/fam) Hiluse vs. ML\n');
unpairedSampleTest(HilusRatio,MlRatio,'verbose',true,'labels',labels);
fprintf('Log2 transformed\n');
%unpairedSampleTest(HilRl2,MlRl2,'verbose',true,'labels',labels);

clear ds allDF outlierThr MlWoD HilWoD signal moving n it labels colors...
    ylims c corrThr f1 outlier_thr pos signaltype smoothFct


%% SUPPLEMENTARY FIGURE 7
%% Preparation: Add FF correlations
% Only matters for the place cell plotting, the other functions consider
% and exclude appropriately

i = 0;
for ds = 1:length(DREADDpv)
    % Add session-stability and trial-trial variability to datasets
    DREADDpv{ds} = reviewer2_corr(DREADDpv{ds});
end

i = 0;
for ds = 1:length(DREADDsom)
    % Add session-stability and trial-trial variability to datasets
    DREADDsom{ds} = reviewer2_corr(DREADDsom{ds});
end

clear i ds

%% A,B) FAM VS NOV COMPARISON DURING BASELINE AND CLOZAPINE
% This version is built for pair-wise comparisons primarily. 
% Cells are selected separately for each comparison, depending on the categories
% selected for comparison.

% SETTING PARAMETERS
plotcats = [1 2; 3 4]; % These are the pairs of categories querried for the violin plots, 
% dim 1 defines number of violins on each subplot
selcats = {[1 2], [3 4]}; % These are the categories used for cell selections 
% (set e.g. to [1 2] to get cells meeting criteria during baseline, [1:6]
% for cells meeting crietria under any condition).
fields = {'transientrate','spatialinfo'}; % This can be any property that is 
% a subfield of 'cells' in the Data structure, e.g 'transientrate','spatial_P',etc.
ACTthr = 1/60;
pSIthr = .05;
pPFthr = 1;
savefigs = false;
combination = 'or';
outlierthr = 3; % Threshold (in interquartile ranges) for tukey outlier removal, set to Inf to disable feature
colors = {[0 0.447 0.741],[0.301 0.745 0.933],...
    [0.85 0.325 0.098],[0.929 0.694 0.125],...
    [.5 .65 .5], [.7 .85 .7]};
xlabels = {'fam-bsl','nov-bsl','fam-clz','nov-clz'};%,'fam-sal','nov-sal'};
dst = .2; % Distance of xticklabels
ylabels = {'Transient rate (Hz)','spatial info (bits/sec)'}; % Needs to match fields
ylims = [-.01 .1; -.1 1.5]; % Plotting range for each parameter. Needs to match fields
sourcenames = {'SupplFig7a','SupplFig7b'};


% EXTRACTING DATA
for f = 1:length(fields)
    % Filling in values for each group of categories.
    for cs = 1:size(plotcats,1)
        PVvalues = metastatistics(DREADDpv,'cell',fields{f},'categories',selcats{cs},'SI',pSIthr,...
            'PF',pPFthr,'active',ACTthr,'combination',combination,'rescategories',plotcats(cs,:),'plotting',false);
        SOMvalues = metastatistics(DREADDsom,'cell',fields{f},'categories',selcats{cs},'SI',pSIthr,...
            'PF',pPFthr,'active',ACTthr,'combination',combination,'rescategories',plotcats(cs,:),'plotting',false);
        
        % Optional: outlier removal here
        PVvalues = tukeyOutlierRemoval(PVvalues,outlierthr);
        SOMvalues = tukeyOutlierRemoval(SOMvalues,outlierthr);

        PVres.(fields{f})(cs,[1 2]) = {PVvalues(:,1),PVvalues(:,2)};
        SOMres.(fields{f})(cs,[1 2]) = {SOMvalues(:,1),SOMvalues(:,2)};
    end
end

% PLOTTING
f1 = figure;
set(f1,'Units','centimeters','defaultAxesFontSize',6,'Position',[4,4,6,10]);

for f = 1:length(fields) 
   % Build 'global' array with all groups for PV and SOM animals
   PVa = catuneven(reshape(PVres.(fields{f})',1,[]),NaN);
   SOMa = catuneven(reshape(SOMres.(fields{f})',1,[]),NaN);

   % Save source data files
   T = array2table(catuneven({PVa,SOMa},NaN),'VariableNames',cat(2,...
       strcat('PV',xlabels),strcat('SOM',xlabels)));
   writetable(T,sourcenames{1,f},'Delimiter','\t','WriteRowNames',false)
   
   % Plot PV data
   %subplot(2,length(fields),f); hold on;
   subplot(length(fields),2,2*f-1);
   stakeplot(PVa,'color',[.3 .3 .3],'LineAlpha',.2,'spread',.45,'jitter',0,'LineWidth',.5);
   distributionPlot(PVres.(fields{f})(:,1),'histOri','left','color',colors(plotcats(:,1)),'widthDiv',[2 1],'showMM',6,'FaceAlpha',.7);
   distributionPlot(PVres.(fields{f})(:,2),'histOri','right','color',colors(plotcats(:,2)),'widthDiv',[2 2],'showMM',6,'FaceAlpha',.7);
   title('PV Cre');
   ylim(ylims(f,:));
   ylabel(ylabels{f});
   xticks(sort([(1:size(plotcats,1))-dst (1:size(plotcats,1))+dst]));
   xticklabels(xlabels(plotcats'));
   xtickangle(45);
   xlim([.4 length(fields)+.6]);
   set(gca,'FontSize',6);
   
   % Plot SOM data
   %subplot(2,length(fields),length(fields)+f);
   subplot(length(fields),2,2*f);
   stakeplot(SOMa,'color',[.3 .3 .3],'LineAlpha',.2,'spread',.45,'jitter',0,'LineWidth',.5);
   distributionPlot(SOMres.(fields{f})(:,1),'histOri','left','color',colors(plotcats(:,1)),'widthDiv',[2 1],'showMM',6,'FaceAlpha',.7);
   distributionPlot(SOMres.(fields{f})(:,2),'histOri','right','color',colors(plotcats(:,2)),'widthDiv',[2 2],'showMM',6,'FaceAlpha',.7);
   title('SOM Cre');
   ylim(ylims(f,:));
   ylabel(ylabels{f});
   xticks(sort([(1:size(plotcats,1))-dst (1:size(plotcats,1))+dst]));
   xticklabels(xlabels(plotcats'));
   xtickangle(45);
   xlim([.4 length(fields)+.6]);
   set(gca,'FontSize',6);
end

% SAVE FIGURE
if savefigs
    pos = get(f1,'Position');
    set(f1,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[pos(3), pos(4)])
    print(f1,sprintf('SInRate.pdf'),'-dpdf','-painters')
end

% STATISTICS
for f = 1:length(fields)
    fprintf('PV Cells %s:\n',fields{f});
    for cs = 1:size(plotcats,1)
        pairedSampleTest(PVres.(fields{f}){cs,1},PVres.(fields{f}){cs,2},...
            'verbose',true,'labels',xlabels(plotcats(cs,:)));
    end
    
    fprintf('SOM Cells %s:\n',fields{f});
    for cs = 1:size(plotcats,1)
        pairedSampleTest(SOMres.(fields{f}){cs,1},SOMres.(fields{f}){cs,2},...
            'verbose',true,'labels',xlabels(plotcats(cs,:)));
    end
end

clear pPFthr pSIthr ACTthr cats colors dst plotcats outlierthr fields...
    ylabels ylims f cs f1 PVvalues SOMvalues pPVs pPVt pSOMs pSOMt...
    PVres SOMres xlabels selcats combination PVa SOMa

%% C) F-N remapping vs. in session stability
labels = {'FF','FN'};
colors = {[.5 .5 .5],[.5 1 .5]};

% PV-Cre DREADD
FFpv = metastatistics(DREADDpv,'cell','sessionstab','active',1/60,'SI',.05,...
    'categories',1,'combination','or','rescategories',1,'plotting',false);
FNpv = metastatistics(DREADDpv,'cell','Pearson_r','active',1/60,'SI',.05,...
    'categories',[1 2],'combination','or','rescategories',[1 2],'plotting',false);
FNpv = FNpv(:,1,2);

figure('position',[200 200 150 220]); hold on;
distributionPlot({FFpv,FNpv},'color',colors,'showMM',6,'distWidth',.9,'addSpread',0,'histOpt',1,'divFactor',3);
scatter(randn(length(FFpv),1)/10+1,FFpv,5,[0 0 0],'.')
scatter(randn(length(FNpv),1)/10+2,FNpv,5,[0 0 0],'.')
title('PV-cre');
ylim([-1.2 1.2]);
ylabel('Spatial Correlation (Pearson''s R)');
xticklabels(labels);
xtickangle(45);

fprintf('PV cells\n')
unpairedSampleTest(FFpv,FNpv,'verbose',true,'labels',labels);

% SOM-Cre DREADD
FFsom = metastatistics(DREADDsom,'cell','sessionstab','active',1/60,'SI',.05,...
    'categories',1,'combination','or','rescategories',1,'plotting',false);
FNsom = metastatistics(DREADDsom,'cell','Pearson_r','active',1/60,'SI',.05,...
    'categories',[1 2],'combination','or','rescategories',[1 2],'plotting',false);
FNsom = FNsom(:,1,2);


figure('position',[200 200 150 220]); hold on;
distributionPlot({FFsom,FNsom},'color',colors,'showMM',6,'distWidth',.9,'addSpread',0,'histOpt',1,'divFactor',3);
scatter(randn(length(FFsom),1)/10+1,FFsom,5,[0 0 0],'.')
scatter(randn(length(FNsom),1)/10+2,FNsom,5,[0 0 0],'.')
title('SOM-cre');
ylim([-1.2 1.2]);
ylabel('Spatial Correlation (Pearson''s R)');
xticklabels(labels);
xtickangle(45);

% Save source data
T = array2table(catuneven({FFpv,FNpv,FFsom,FNsom},NaN),'VariableNames',...
    {'FF corr PV','FN corr PV','FF corr SOM','FN corr SOM'});
writetable(T,'SupplFig7c','Delimiter','\t','WriteRowNames',false)

fprintf('SOM cells\n')
unpairedSampleTest(FFsom,FNsom,'verbose',true,'labels',labels);

clear FFpv FNpv FFsom FNsom colors labels

%% D) Population vector analysis
datasets = {'DREADDpv','DREADDsom'};
cats = [1 2; 3 4];
colors = [0 0.447 0.741 ; 0.85 0.325 0.098];%0.301 0.745 0.933];% ;  ; 0.929 0.694 0.125];
f1 = figure('position',[400 400 150 250]); hold on;
sourcedata = {};

for nmds = 1:length(datasets)
    clear FNcorr
    eval(sprintf('mds = %s;',datasets{nmds}));
    for ds = 1:length(mds)
        for c = [1,2]
            selected = findcells(mds{ds},cats(c,:),1/60,.05,...
                1,'or','');
            if length(selected)>4
                thisDS = PVcorrelation(mds{ds},selected);
                FNcorr(ds,c) = thisDS.metadata.PVcorr(cats(c,1),cats(c,2));
            else
                FNcorr(ds,c) = NaN;
            end
        end
    end
    
    % Plotting
    xv = [(nmds-1)*2 (nmds-1)*2+1];
    bA = bar(xv,nanmean(FNcorr,1),'FaceColor','flat');
    bA.CData = colors;
    bA.EdgeColor= [.5 .5 .5];
    bA.LineWidth=.5;
    
    plot(xv,FNcorr,'-','color',[.2 .2 .2],'LineWidth',.5);
    ylabel('Population Vector Corr'); 
    
    %title(datasets{nmds});
    
    % Statistics
    fprintf('%s statistics:\n',datasets{nmds});
    pairedSampleTest(FNcorr(:,1),FNcorr(:,2),'verbose',true);
    
    sourcedata{end+1} = FNcorr(:,1);
    sourcedata{end+1} = FNcorr(:,2);
end

% Save source data
T = array2table(catuneven(sourcedata,NaN),'VariableNames',...
    {'PV bsl','PV clz','SOM bsl','SOM clz'});
writetable(T,'SupplFig7d','Delimiter','\t','WriteRowNames',false)

xticks(0:3); xtickangle(45); xticklabels({'Bsl','Clz','Bsl','Clz'});
ylim([-.01 .7]);

% SAVE FIGURE
pos = get(f1,'Position');
set(f1,'PaperPositionMode','Auto','PaperUnits','points','PaperSize',[pos(3), pos(4)])
print(f1,'PopVectCorr.pdf','-dpdf','-painters')


clear nmds ds mds datasets FNcorr thisDS cats c selected


%% SUPPLEMENTARY FIGURE 8  -- CA1 SOM and CA1 PV DREADD MICE
%% B,C) Rates and SI CA1 PV
% This version is built for pair-wise comparisons primarily.

% SETTING PARAMETERS
plotcats = [1 3; 2 4]; % These are the pairs of categories querried for the violin plots, 
% dim 1 defines number of violins on each subplot
selcats = {[1 3], [2 4]}; % These are the categories used for cell selections 
% (set e.g. to [1 2] to get cells meeting criteria during baseline, [1:6]
% for cells meeting crietria under any condition).
fields = {'transientrate','spatialinfo'}; % This can be any property that is 
% a subfield of 'cells' in the Data structure, e.g 'transientrate','spatial_P',etc.
ACTthr = 1/60;
pSIthr = .05;
pPFthr = 1;
combination = 'or';
outlierthr = 3; % Threshold (in interquartile ranges) for tukey outlier removal, set to Inf to disable feature
colors = {[0 0.447 0.741],[0.301 0.745 0.933],...
    [0.85 0.325 0.098],[0.929 0.694 0.125],...
    [.5 .65 .5], [.7 .85 .7]};
xlabels = {'fam-bsl','nov-bsl','fam-clz','nov-clz'};
dst = .2; % Distance of xticklabels
ylabels = {'Transient rate (Hz)','spatial info (bits/sec)'}; % Needs to match fields
ylims = [-.01 .15; -.4 7]; % Plotting range for each parameter. Needs to match fields
sourcenames = {'SupplFig8b','SupplFig8c'};


% EXTRACTING DATA
for f = 1:length(fields)
    % Filling in values for each group of categories.
    for cs = 1:size(plotcats,1)
        CA1PVvalues = metastatistics(DREADDca1PV,'cell',fields{f},'categories',selcats{cs},'SI',pSIthr,...
            'PF',pPFthr,'active',ACTthr,'combination',combination,'rescategories',plotcats(cs,:),'plotting',false);
        
        % Optional: outlier removal here
        CA1PVvalues = tukeyOutlierRemoval(CA1PVvalues,outlierthr);

        CA1PVres.(fields{f})(cs,[1 2]) = {CA1PVvalues(:,1),CA1PVvalues(:,2)};
    end
end

% PLOTTING
f1 = figure;
set(f1,'Units','centimeters','defaultAxesFontSize',8,'Position',[4,4,7,5]);

for f = 1:length(fields) 
   % Build 'global' array with all groups for PV and SOM animals
   CA1PVa = catuneven(reshape(CA1PVres.(fields{f})',1,[]),NaN);
   
   % Save source data files
   T = array2table(CA1PVa,'VariableNames',xlabels(plotcats'));
   writetable(T,sourcenames{f},'Delimiter','\t','WriteRowNames',false)
   
   subplot(1,length(fields),f);
   stakeplot(CA1PVa,'color',[.2 .2 .2],'LineAlpha',.1,'spread',.45,'jitter',0,'LineWidth',.5);
   distributionPlot(CA1PVres.(fields{f})(:,1),'histOri','left','color',colors(plotcats(:,1)),'widthDiv',[2 1],'showMM',6,'FaceAlpha',.7);
   distributionPlot(CA1PVres.(fields{f})(:,2),'histOri','right','color',colors(plotcats(:,2)),'widthDiv',[2 2],'showMM',6,'FaceAlpha',.7);
   title('CA1 PV Cre');
   ylim(ylims(f,:));
   ylabel(ylabels{f});
   xticks(sort([(1:size(plotcats,1))-dst (1:size(plotcats,1))+dst]));
   xticklabels(xlabels(plotcats'));
   xtickangle(45);
   xlim([.4 length(fields)+.6]);
end

% SAVE FIGURE
pos = get(f1,'Position');
set(f1,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[pos(3), pos(4)])
print(f1,sprintf('PV_SInRate.pdf'),'-dpdf','-painters')

% STATISTICS
for f = 1:length(fields)
    fprintf('PV Cells %s:\n',fields{f});
    for cs = 1:size(plotcats,1)
        pairedSampleTest(CA1PVres.(fields{f}){cs,1},CA1PVres.(fields{f}){cs,2},...
            'verbose',true,'labels',xlabels(plotcats(cs,:)));
    end
end

clear pPFthr pSIthr ACTthr cats colors dst plotcats outlierthr fields...
    ylabels ylims f cs f1 PVvalues SOMvalues pPVs pPVt pSOMs pSOMt...
    PVres SOMres xlabels selcats combination PVa SOMa

%% D) FN Remapping BSL VS CLZ CA1 PV
labels = {'bsl FN','clz FN'};
colors = {[0 0.447 0.741],[0.85 0.325 0.098]};

[res,~] = metacorr_pairwise(DREADDca1PV,'categories',1:4,'active',1/60,...
    'SI',.05,'combination','or','plotresults',false);

FNbsl = squeeze(res(1,2,:));
FNclz = squeeze(res(3,4,:));

f1 = figure; hold on;
set(f1,'Units','centimeters','defaultAxesFontSize',8,'Position',[4,4,3.5,5]);

distributionPlot({FNbsl,FNclz},'color',colors,'showMM',6,'distWidth',.9,...
    'addSpread',0,'histOpt',1,'divFactor',3,'FaceAlpha',.5);
scatter(randn(length(FNbsl),1)/10+1,FNbsl,8,max(colors{1}-.2,0),'.')
scatter(randn(length(FNclz),1)/10+2,FNclz,8,max(colors{2}-.2,0),'.')
ylim([-.7 1.2]);
ylabel('Spatial Correlation (R)');
xticklabels(labels);
xtickangle(45);
set(gca,'FontSize',6);

% Save source data files
T = array2table(catuneven({FNbsl,FNclz},NaN),'VariableNames',...
    {'FN corr bsl','FN corr clz'});
writetable(T,'SupplFig8d','Delimiter','\t','WriteRowNames',false)

% SAVE FIGURE
pos = get(f1,'Position');
set(f1,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[pos(3), pos(4)])
print(f1,sprintf('PV_FN_bslVSclz.pdf'),'-dpdf','-painters','-bestfit')

fprintf('Comparing PV cells with place fields during bsl, or clz respectively:\n')
unpairedSampleTest(FNbsl,FNclz,'verbose',true,'labels',labels);

clear res labels FNbsl FNclz

%% E) Place maps CA1 PV
plot_placemaps_meta(DREADDca1PV,1,'categories',[1 2],'SI',.05,'active',1/60,'combination','or','clims',[0 .5],'figscale',[.5 .5]);
plot_placemaps_meta(DREADDca1PV,3,'categories',[3 4],'SI',.05,'active',1/60,'combination','or','clims',[0 .5],'figscale',[.5 .5]);

%% G,H) Rates and SI CA1 SOM
% This version is built for pair-wise comparisons primarily.

% SETTING PARAMETERS
plotcats = [1 3; 2 4]; % These are the pairs of categories querried for the violin plots, 
% dim 1 defines number of violins on each subplot
selcats = {[1 3], [2 4]}; % These are the categories used for cell selections 
% (set e.g. to [1 2] to get cells meeting criteria during baseline, [1:6]
% for cells meeting crietria under any condition).
fields = {'transientrate','spatialinfo'}; % This can be any property that is 
% a subfield of 'cells' in the Data structure, e.g 'transientrate','spatial_P',etc.
ACTthr = 1/60;
pSIthr = .05;
pPFthr = 1;
combination = 'or';
outlierthr = 3; % Threshold (in interquartile ranges) for tukey outlier removal, set to Inf to disable feature
colors = {[0 0.447 0.741],[0.301 0.745 0.933],...
    [0.85 0.325 0.098],[0.929 0.694 0.125],...
    [.5 .65 .5], [.7 .85 .7]};
xlabels = {'fam-bsl','nov-bsl','fam-clz','nov-clz'};
dst = .2; % Distance of xticklabels
ylabels = {'Transient rate (Hz)','spatial info (bits/sec)'}; % Needs to match fields
ylims = [-.01 .3; -.08 3]; % Plotting range for each parameter. Needs to match fields
sourcenames = {'SupplFig8g','SupplFig8h'};


% EXTRACTING DATA
for f = 1:length(fields)
    % Filling in values for each group of categories.
    for cs = 1:size(plotcats,1)
        SOMvalues = metastatistics(DREADDca1SOM,'cell',fields{f},'categories',selcats{cs},'SI',pSIthr,...
            'PF',pPFthr,'active',ACTthr,'combination',combination,'rescategories',plotcats(cs,:),'plotting',false);
        
        % Optional: outlier removal here
        SOMvalues = tukeyOutlierRemoval(SOMvalues,outlierthr);

        SOMres.(fields{f})(cs,[1 2]) = {SOMvalues(:,1),SOMvalues(:,2)};
    end
end

% PLOTTING
f1 = figure;
set(f1,'Units','centimeters','defaultAxesFontSize',8,'Position',[4,4,7,5]);

for f = 1:length(fields) 
   % Build 'global' array with all groups for PV and SOM animals
   SOMa = catuneven(reshape(SOMres.(fields{f})',1,[]),NaN);

   % Save source data files
   T = array2table(SOMa,'VariableNames',xlabels(plotcats'));
   writetable(T,sourcenames{f},'Delimiter','\t','WriteRowNames',false)
   
   subplot(1,length(fields),f);
   stakeplot(SOMa,'color',[.2 .2 .2],'LineAlpha',.1,'spread',.45,'jitter',0,'LineWidth',.5);
   distributionPlot(SOMres.(fields{f})(:,1),'histOri','left','color',colors(plotcats(:,1)),'widthDiv',[2 1],'showMM',6,'FaceAlpha',.7);
   distributionPlot(SOMres.(fields{f})(:,2),'histOri','right','color',colors(plotcats(:,2)),'widthDiv',[2 2],'showMM',6,'FaceAlpha',.7);
   title('SOM Cre');
   ylim(ylims(f,:));
   ylabel(ylabels{f});
   xticks(sort([(1:size(plotcats,1))-dst (1:size(plotcats,1))+dst]));
   xticklabels(xlabels(plotcats'));
   xtickangle(45);
   xlim([.4 length(fields)+.6]);
end

% SAVE FIGURE
pos = get(f1,'Position');
set(f1,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[pos(3), pos(4)])
print(f1,sprintf('SOM_CA1_SInRate.pdf'),'-dpdf','-painters')

% STATISTICS
for f = 1:length(fields)
    fprintf('SOM Cells %s:\n',fields{f});
    for cs = 1:size(plotcats,1)
        pairedSampleTest(SOMres.(fields{f}){cs,1},SOMres.(fields{f}){cs,2},...
            'verbose',true,'labels',xlabels(plotcats(cs,:)));
    end
end

clear pPFthr pSIthr ACTthr cats colors dst plotcats outlierthr fields...
    ylabels ylims f cs f1 PVvalues SOMvalues pPVs pPVt pSOMs pSOMt...
    PVres SOMres xlabels selcats combination PVa SOMa

%% I) FN Remapping BSL VS CLZ CA1 SOM
labels = {'bsl FN','clz FN'};
colors = {[0 0.447 0.741],[0.85 0.325 0.098]};

[res,~] = metacorr_pairwise(DREADDca1SOM,'categories',1:4,'active',1/60,...
    'SI',.05,'combination','or','plotresults',false);

FNbsl = squeeze(res(1,2,:));
FNclz = squeeze(res(3,4,:));

f1 = figure; hold on;
set(f1,'Units','centimeters','defaultAxesFontSize',8,'Position',[4,4,3.5,5]);

distributionPlot({FNbsl,FNclz},'color',colors,'showMM',6,'distWidth',.9,...
    'addSpread',0,'histOpt',1,'divFactor',3,'FaceAlpha',.5);
scatter(randn(length(FNbsl),1)/10+1,FNbsl,8,max(colors{1}-.2,0),'.')
scatter(randn(length(FNclz),1)/10+2,FNclz,8,max(colors{2}-.2,0),'.')
ylim([-.7 1.2]);
ylabel('Spatial Correlation (R)');
xticklabels(labels);
xtickangle(45);
set(gca,'FontSize',6);

% Save source data files
T = array2table(catuneven({FNbsl,FNclz},NaN),'VariableNames',...
    {'FN corr bsl','FN corr clz'});
writetable(T,'SupplFig8i','Delimiter','\t','WriteRowNames',false)

% SAVE FIGURE
pos = get(f1,'Position');
set(f1,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[pos(3), pos(4)])
print(f1,sprintf('SOM_FN_bslVSclz.pdf'),'-dpdf','-painters','-bestfit')

fprintf('Comparing SOM cells with place fields during bsl, or clz respectively:\n')
unpairedSampleTest(FNbsl,FNclz,'verbose',true,'labels',labels);

clear res labels FNbsl FNclz

%% J) Place maps CA1 SOM
plot_placemaps_meta(DREADDca1SOM,1,'categories',[1 2],'SI',.05,'active',1/60,'combination','or','clims',[0 .5],'figscale',[.5 .5]);
plot_placemaps_meta(DREADDca1SOM,3,'categories',[3 4],'SI',.05,'active',1/60,'combination','or','clims',[0 .5],'figscale',[.5 .5]);

%% K,L) 3-way ANOVA on Pre-post clz ratios (Region, SOM/PV, FamNov)
% SETTING PARAMETERS
plotcats = [1 3; 2 4]; % These are the pairs of categories querried for the violin plots, 
% dim 1 defines number of violins on each subplot
selcats = {[1 3], [2 4]}; % These are the categories used for cell selections 
% (set e.g. to [1 2] to get cells meeting criteria during baseline, [1:6]
% for cells meeting crietria under any condition).
fields = {'transientrate','spatialinfo'}; % This can be any property that is 
datasets = {'DREADDpv','DREADDsom','DREADDca1PV','DREADDca1SOM'};
grouplabels = {'DG PV Fam','DG PV Nov','DG SOM Fam','DG SOM Nov',...
        'CA1 PV Fam','CA1 PV Nov','CA1 SOM Fam','CA1 SOM Nov'};
% a subfield of 'cells' in the Data structure, e.g 'transientrate','spatial_P',etc.
ACTthr = 1/60;
pSIthr = .05;
pPFthr = 1;
savefigs = false;
combination = 'or';
outlierthr = 3; % Threshold (in interquartile ranges) for tukey outlier removal, set to Inf to disable feature
colors = repmat({[0.85 0.325 0.098],[0.929 0.694 0.125]},1,4);%[0 0.447 0.741],[0.301 0.745 0.933],...
    %[.5 .65 .5], [.7 .85 .7]};
sourcenames = {'SupplFig8k','SupplFig8l'};

res = {[],[]};

% EXTRACTING DATA
for f = 1:length(fields)
    for ds = 1:length(datasets)
        
        eval(sprintf('metads = %s;',datasets{ds}));
        
        for cs = 1:size(plotcats,1)
            values = metastatistics(metads,'cell',fields{f},'categories',selcats{cs},'SI',pSIthr,...
                'PF',pPFthr,'active',ACTthr,'combination',combination,'rescategories',plotcats(cs,:),'plotting',false);
            
            ratios = values(:,2)./values(:,1);
            
            % Optional: outlier removal here
            ratios = tukeyOutlierRemoval(ratios,outlierthr);
            
            resCell{f,cs,ds} = ratios;
            
            % ANOVA labels (column 2 = PV/SOM, 3 = Fam/Nov
            ratios(:,2) = mod(ds+1,2); % Column 2: PV = 0, SOM = 1
            ratios(:,3) = cs-1; % Column 3: 0 = fam, 1 = nov
            ratios(:,4) = double(ds>2); % 0 = DG, 1 = CA1;
            
            res{f} = cat(1,res{f},ratios);
        end
    end
    
    % PLOT RATIOS
    f1 = figure; hold on;
    plotvals = real(log2(res{f}(:,1)));
    [~,isoutlier] = tukeyOutlierRemoval(plotvals,3);
    
    set(f1,'Units','centimeters','defaultAxesFontSize',10,'Position',[4,4,8.6,6.6]);
    
    % Group index = PVfam(0), PVnov(1), SOMfam(0+2), SOMnov(1+2);
    distributionPlot(plotvals(~isoutlier),'groups',...
        res{f}(~isoutlier,3) + res{f}(~isoutlier,2)*2 + res{f}(~isoutlier,4)*4,...
        'color',colors,'showMM',6,'distWidth',.9,...
        'addSpread',0,'histOpt',1,'divFactor',3,'FaceAlpha',.5);
    
    plotSpread(plotvals(~isoutlier),'distributionIdx',...
        res{f}(~isoutlier,3) + res{f}(~isoutlier,2)*2 + res{f}(~isoutlier,4)*4,...
        'spreadFcn',{'lin',10},'spreadWidth',.9,...
        'distributionColors',repmat({[max(colors{1}-.2,0) .01],[max(colors{2}-.2,0) .01]},1,4),...
        'distributionMarkerSize',1,'binWidth',.05);
        
    ylabel(sprintf('Log2 %s ratio clz/bsl',fields{f}));
    xticks(0:7);
    xticklabels(grouplabels);
    xtickangle(45);
    
    pos = get(f1,'Position');
    set(f1,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[pos(3), pos(4)])
    print(f1,sprintf('%s ratios 3way.pdf',fields{f}),'-dpdf','-painters')    

    % Save source data files
    T = array2table(catuneven(cat(3,resCell(f,1,:),resCell(f,2,:)),NaN),'VariableNames',...
        cat(2,strcat({'DG PV','DG SOM','CA1 PV','CA1 SOM'},' FAM'),...
        strcat({'DG PV','DG SOM','CA1 PV','CA1 SOM'},' NOV')));
    writetable(T,sourcenames{f},'Delimiter','\t','WriteRowNames',false)
    
    % 3-WAY ANOVA
    [~,tbl,stats] = anovan(res{f}(:,1),{res{f}(:,2),res{f}(:,3),res{f}(:,4)},...
        "Model","interaction","Varnames",["SOM","Nov","CA1"]);
    fprintf('3-way ANOVA with interactions for %s in SOM-PV, FAM-NOV:\n',fields{f});
    fprintf('%s\t%s\t%s\t%s\t%s\t%s\t%s\n',tbl{1,1},tbl{1,2},tbl{1,3},tbl{1,4},tbl{1,5},tbl{1,6},tbl{1,7})
    for n = 2:size(tbl,1)
        fprintf('%s\t%.2e\t%i\t%.2g\t%.2f\t%.2g\t%.4g\n',tbl{n,1},tbl{n,2},tbl{n,3},tbl{n,4},tbl{n,5},tbl{n,6},tbl{n,7})
    end
    
    figure;
    [compc,~,~,gnames] = multcompare(stats,"Dimension",[1 2 3]);
    fprintf('\nDunn\''s post-hoc test\n')
    fprintf('%s\t%s\t%s\n','Group1','Group2','p')
    for n = 1:size(compc,1)
        fprintf('%s\t%s\t%.2e\n',gnames{compc(n,1)},gnames{compc(n,2)},compc(n,6))
    end
end

clear pPFthr pSIthr ACTthr cats colors dst plotcats outlierthr fields...
    ylabels ylims f cs f1 values ratios...
    xlabels selcats combination PVa SOMa metads

