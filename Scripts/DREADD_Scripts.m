%% Scripts for DREADD and axon imaging data by Aurore Cazala
% Load 'DataPV Aurore Cazala.mat' and 'Data_Som Aurore Cazala (1).mat',
% rename 'Data' to 'DREADDpv' and 'DREADDsom', respectively.
% 
% Infos:
% - Categories 1,2 = Fam/Nov baseline, 3,4 = clozapine, [5,6 = saline]

%% OBSOLETE === Get place cell numbers per experiments in table
for ds = 1:length(DREADDpv)
    for n = 1:length(DREADDpv{ds}.cells)
        isplc(n) = any(DREADDpv{ds}.cells{n}.Placefield_P(1:4)<0.05);
    end
    PVnPlc(ds) = sum(isplc);
    %isplc = DREADDpv{ds}.metadata.Placefield_P(:,1:4)<0.05;
    %PVnPlc(ds) = sum(any(isplc,2));
end
sum(PVnPlc)

for ds = 1:length(DREADDsom)
    for n = 1:length(DREADDsom{ds}.cells)
        isplc(n) = any(DREADDsom{ds}.cells{n}.Placefield_P(1:4)<0.05);
    end
    SOMnPlc(ds) = sum(isplc);
    %isplc = DREADDsom{ds}.metadata.Placefield_P(:,1:4)<0.05;
    %SOMnPlc(ds) = sum(any(isplc,2));
end
sum(SOMnPlc)

clear isplc ds n

%% OBSOLETE === Reverse-engineer place cells from rates / SI fields in AC data
% - The maps on the figures can be re-created by using 'SI' < 0.05 and
% 'active' < 1/60 for the category to be plotted! Limit caxis([0 .2]);

%% Recalculate place fields
for ds = 1:length(DREADDpv)
    DREADDpv{ds} = placefields(DREADDpv{ds});
end

for ds = 1:length(DREADDsom)
    DREADDsom{ds} = placefields(DREADDsom{ds});
end

plot_placemaps_meta(DREADDsom,1,'categories',1:4,'PF',0.05,'combination','or');
plot_placemaps_meta(DREADDpv,1,'categories',1:4,'PF',0.05,'combination','or');

clear ds

%% Plot rates and SI (violin plots) - flexible number of variables, categories and ability to querry nov/fam cells independently
% Setting parameters
cats = [1 2; 3 4]; % The categories between semicola are querried separately, set to [1 3; 2 4] to get fam/nov place cells separately, set 1:4 for all place cells
plotcats = [1 2; 3 4]; % Define the order for plotting
fields = {'transientrate','spatialinfo'}; % This can be any property that is a subfield of 'cells' in the Data structure, e.g 'transientrate','spatial_P',etc.
ACTthr = 1/60;
pSIthr = .05;
pPFthr = 1;
outlierthr = 3; % Threshold (in interquartile ranges) for tukey outlier removal, set to Inf to disable feature
colors = {[0 0.447 0.741],[0.301 0.745 0.933],...
    [0.85 0.325 0.098],[0.929 0.694 0.125]};
xlabels = {'fam-bsl','nov-bsl','fam-clz','nov-clz'};
dst = .2; % Distance of xticklabels
ylabels = {'Transient rate (Hz)','spatial info (bits/sec)'}; % Needs to match fields
ylims = [-.01 .1; -.1 1.5]; % Plotting range for each parameter. Needs to match fields


% EXTRACTING DATA
for f = 1:length(fields)
    % Initializing resulttable as large NaN array, 1st dimension needs to be well above any expected cell number
    PVres.(fields{f}) = NaN(1e4,numel(cats)); 
    SOMres.(fields{f}) = NaN(1e4,numel(cats));
    
    % Filling in values for each group of categories.
    for cs = 1:size(cats,1)
        PVvalues = metastatistics(DREADDpv,'cell',fields{f},'categories',cats(cs,:),'SI',pSIthr,...
            'PF',pPFthr,'active',ACTthr,'combination','or','rescategories',cats(cs,:),'plotting',false);
        SOMvalues = metastatistics(DREADDsom,'cell',fields{f},'categories',cats(cs,:),'SI',pSIthr,...
            'PF',pPFthr,'active',ACTthr,'combination','or','rescategories',cats(cs,:),'plotting',false);
        
        % Optional: outlier removal here
        PVvalues = tukeyOutlierRemoval(PVvalues,outlierthr);
        SOMvalues = tukeyOutlierRemoval(SOMvalues,outlierthr);
        
        PVres.(fields{f})(1:size(PVvalues,1),cats(cs,:)) = PVvalues;
        SOMres.(fields{f})(1:size(SOMvalues,1),cats(cs,:)) = SOMvalues;
    end
    
    % Removing excess rows at the bottom of the array.
    PVres.(fields{f}) = PVres.(fields{f})(~all(isnan(PVres.(fields{f})),2),:);
    SOMres.(fields{f}) = SOMres.(fields{f})(~all(isnan(SOMres.(fields{f})),2),:);
end

% PLOTTING
f1 = figure;
for f = 1:length(fields) 
   % Plot PV data
   subplot(2,length(fields),f);
   distributionPlot(PVres.(fields{f})(:,plotcats(:,1)),'histOri','left','color',colors(plotcats(:,1)),'widthDiv',[2 1],'showMM',6);
   distributionPlot(PVres.(fields{f})(:,plotcats(:,2)),'histOri','right','color',colors(plotcats(:,2)),'widthDiv',[2 2],'showMM',6);
   title('PV Cre');
   ylim(ylims(f,:));
   ylabel(ylabels{f});
   xticks([1-dst 1+dst 2-dst 2+dst]);
   xticklabels(xlabels(plotcats'));
   xtickangle(45);
   
   subplot(2,length(fields),length(fields)+f);
   distributionPlot(SOMres.(fields{f})(:,plotcats(:,1)),'histOri','left','color',colors(plotcats(:,1)),'widthDiv',[2 1],'showMM',6);
   distributionPlot(SOMres.(fields{f})(:,plotcats(:,2)),'histOri','right','color',colors(plotcats(:,2)),'widthDiv',[2 2],'showMM',6);
   title('SOM Cre');
   ylim(ylims(f,:));
   ylabel(ylabels{f});
   xticks([1-dst 1+dst 2-dst 2+dst]);
   xticklabels(xlabels(plotcats'));
   xtickangle(45);
end

% STATISTICS
for f = 1:length(fields)
    fprintf('PV Cells %s:\n',fields{f});
    pairedSampleTest(PVres.(fields{f})(:,cats(1,1)),PVres.(fields{f})(:,cats(1,2)),...
        'verbose',true,'labels',xlabels(cats(1,:)));
    pairedSampleTest(PVres.(fields{f})(:,cats(2,1)),PVres.(fields{f})(:,cats(2,2)),...
        'verbose',true,'labels',xlabels(cats(2,:)));
    
    fprintf('SOM Cells %s:\n',fields{f});
    pairedSampleTest(SOMres.(fields{f})(:,cats(1,1)),SOMres.(fields{f})(:,cats(1,2)),...
        'verbose',true,'labels',xlabels(cats(1,:)));
    pairedSampleTest(SOMres.(fields{f})(:,cats(2,1)),SOMres.(fields{f})(:,cats(2,2)),...
        'verbose',true,'labels',xlabels(cats(2,:)));
end

% kruskalwallisNdunns(PVres.spatialinfo,xlabels); -- Don't (not the same
% groups)
% [h,p] = ttest(PVres.spatialinfo(:,1),PVres.spatialinfo(:,3))
% [p,h] = signrank(SOMres.transientrate(:,1),SOMres.transientrate(:,3))

clear pPFthr pSIthr ACTthr cats colors dst plotcats outlierthr fields...
    ylabels ylims f cs f1 PVvalues SOMvalues pPVs pPVt pSOMs pSOMt...
    PVres SOMres xlabels


%% Plot rates and SI (violin plots) - simple initial version, selects only one group of cells
% Setting parameters
cats = [2 4]; % Categories for place cell / parameter checking. Will always plot 4 categories
ACTthr = 1/60;
pSIthr = .05;
pPFthr = 1;
colors = {[0 0.447 0.741],[0.301 0.745 0.933],[0.85 0.325 0.098],[0.929 0.694 0.125]};
labels = {'fam-bsl','nov-bsl','fam-clz','nov-clz'};

% Extracting Data
PVrates = metastatistics(DREADDpv,'cell','transientrate','categories',cats,'SI',pSIthr,...
    'PF',pPFthr,'active',ACTthr,'combination','or','rescategories',1:4,'plotting',false);
SOMrates = metastatistics(DREADDsom,'cell','transientrate','categories',cats,'SI',pSIthr,...
    'PF',pPFthr,'active',ACTthr,'combination','or','rescategories',1:4,'plotting',false);

PVsi = metastatistics(DREADDpv,'cell','spatialinfo','categories',cats,'SI',pSIthr,...
    'PF',pPFthr,'active',ACTthr,'combination','or','rescategories',1:4,'plotting',false);
SOMsi = metastatistics(DREADDsom,'cell','spatialinfo','categories',cats,'SI',pSIthr,...
    'PF',pPFthr,'active',ACTthr,'combination','or','rescategories',1:4,'plotting',false);

% Plotting
f1 = figure;
subplot(2,2,1)
distributionPlot(PVrates(:,[1 2]),'histOri','left','color',colors([1 2]),'widthDiv',[2 1],'showMM',6);
distributionPlot(PVrates(:,[3 4]),'histOri','right','color',colors([3 4]),'widthDiv',[2 2],'showMM',6);
title('DREADDs PV-Cre')
ylabel('transient rates (Hz)');
ylim([0 .1]);

subplot(2,2,2)
distributionPlot(SOMrates(:,[1 2]),'histOri','left','color',colors([1 2]),'widthDiv',[2 1],'showMM',6);
distributionPlot(SOMrates(:,[3 4]),'histOri','right','color',colors([3 4]),'widthDiv',[2 2],'showMM',6);
title('DREADDs SOM-Cre')
ylabel('transient rates (Hz)');
ylim([0 .1]);

subplot(2,2,3)
distributionPlot(PVsi(:,[1 2]),'histOri','left','color',colors([1 2]),'widthDiv',[2 1],'showMM',6);
distributionPlot(PVsi(:,[3 4]),'histOri','right','color',colors([3 4]),'widthDiv',[2 2],'showMM',6);
title('DREADDs PV-Cre')
ylabel('Spatial info');
ylim([0 1]);

subplot(2,2,4)
distributionPlot(SOMsi(:,[1 2]),'histOri','left','color',colors([1 2]),'widthDiv',[2 1],'showMM',6);
distributionPlot(SOMsi(:,[3 4]),'histOri','right','color',colors([3 4]),'widthDiv',[2 2],'showMM',6);
title('DREADDs SOM-Cre')
ylabel('Spatial info');
ylim([0 1]);

clear pPFthr pSIthr ACTthr cats colors
% clear PVsi SOMsi PVrates SOMrates
% kruskalwallisNdunns(PVrates,labels);

%% F-statistic to assess for unequal variances (pre-post clozapine)
% F = ratio mean variance clz / mean variance bsl - larger variance goes to
% numerator!; Degrees of freedom = sample size (N) -1, respectively. Take
% p*2 for 2-tailed test!
%
% p = (1-fpdf(F,N1-1,N2-1))*2;
cats = [1 3];

PVfplc = metastatistics(DREADDpv,'cell','transientrate','categories',cats,'SI',.05,...
    'PF',1,'active',1/60,'combination','or','rescategories',cats,'plotting',false);
SOMfplc = metastatistics(DREADDsom,'cell','transientrate','categories',cats,'SI',.05,...
    'PF',1,'active',1/60,'combination','or','rescategories',cats,'plotting',false);

Fpv = var(PVfplc(:,2))/var(PVfplc(:,1)); if Fpv<1; Fpv = 1/Fpv; end
Fsom = var(SOMfplc(:,2))/var(SOMfplc(:,1)); if Fsom<1; Fsom = 1/Fsom; end

pPV = (1-fcdf(Fpv,size(PVfplc,1)-1,size(PVfplc,1)-1))*2;
pSOM = (1-fcdf(Fsom,size(SOMfplc,1)-1,size(SOMfplc,1)-1))*2;

fprintf('Comparison categories %i vs. %i\n',cats(1),cats(2));
fprintf('PV: F = %.2f, p = %.2e, %i degrees of freedom\n',Fpv,pPV,size(PVfplc,1)-1)
fprintf('SOM: F = %.2f, p = %.2e, %i degrees of freedom\n',Fsom,pSOM,size(SOMfplc,1)-1)

clear PVfplc SOMfplc cats

%% Compare fractions of place cells between fam / nov and bsl / clz
cats = [2 4]; % Categories to compare (matters only for the stats section)
labels = {'fam bsl','nov bsl','fam clz','nov clz'};
PltOrd = [1 3 2 4]; % Order of categories in bars
colors = [0 0.447 0.741 ; 0.301 0.745 0.933 ; 0.85 0.325 0.098 ; 0.929 0.694 0.125];

% GET PV DATA
PVnplc = metastatistics(DREADDpv,'dataset','ncells','categories',1:4,...
    'rescategories',1:4,'active',1/60,'SI',.05,'plotting',false);
PVnact = metastatistics(DREADDpv,'dataset','ncells','categories',1:4,...
    'rescategories',1:4,'active',1/60,'SI',1,'plotting',false);
PVntot = metastatistics(DREADDpv,'dataset','ncells','categories',1:4,...
    'rescategories',1:4,'plotting',false);

PVfplc = PVnplc ./ PVntot .* 100;
PVfact = PVnact ./ PVntot .* 100;

PVcontMtrx = [sum(PVnplc(:,cats(1))) sum(PVntot(:,cats(1)))-sum(PVnplc(:,cats(1)));...
    sum(PVnplc(:,cats(2))) sum(PVntot(:,cats(2)))-sum(PVnplc(:,cats(2)))];


% GET SOM DATA
SOMnplc = metastatistics(DREADDsom,'dataset','ncells','categories',1:4,...
    'rescategories',1:4,'active',1/60,'SI',.05,'plotting',false);
SOMnact = metastatistics(DREADDsom,'dataset','ncells','categories',1:4,...
    'rescategories',1:4,'active',1/60,'SI',1,'plotting',false);
SOMntot = metastatistics(DREADDsom,'dataset','ncells','categories',1:4,...
    'rescategories',1:4,'plotting',false);

SOMfplc = SOMnplc ./ SOMntot .* 100;
SOMfact = SOMnact ./ SOMntot .* 100;

SOMcontMtrx = [sum(SOMnplc(:,cats(1))) sum(SOMntot(:,cats(1)))-sum(SOMnplc(:,cats(1)));...
    sum(SOMnplc(:,cats(2))) sum(SOMntot(:,cats(2)))-sum(SOMnplc(:,cats(2)))];

% PLOTTING
figure
subplot(1,2,1); hold on
bA=bar(sum(PVnact(:,PltOrd))./sum(PVntot).*100,'FaceColor','flat'); % CAVE: This is DIFFERENT from the average of the fractions!
bP=bar(sum(PVnplc(:,PltOrd))./sum(PVntot).*100,'FaceColor','flat');
bA.CData = colors(PltOrd,:);
bP.CData = colors(PltOrd,:);
bP.LineStyle='--';bP.LineWidth=1; bA.EdgeColor= [.5 .5 .5];
bA.LineWidth=1; alpha(.6); ylabel('cell percent'); xticks(1:4); xtickangle(45);

plot([1 2],PVfplc(:,PltOrd([1 2]))','-','color',[.2 .2 .2],'LineWidth',1);
plot([3 4],PVfplc(:,PltOrd([3 4]))','-','color',[.2 .2 .2],'LineWidth',1);

title('DREADD PV');
set(gca,'xticklabel',labels(PltOrd),'Fontsize',20,'fontname','arial','Fontweight','normal',...
    'TickDir','out','lineWidth',2,'box','off');

subplot(1,2,2); hold on
bA=bar(sum(SOMnact(:,PltOrd))./sum(SOMntot).*100,'FaceColor','flat');
bP=bar(sum(SOMnplc(:,PltOrd))./sum(SOMntot).*100,'FaceColor','flat');
bA.CData = colors(PltOrd,:);
bP.CData = colors(PltOrd,:);
bP.LineStyle='--';bP.LineWidth=1; bA.EdgeColor= [.5 .5 .5];
bA.LineWidth=1; alpha(.6); ylabel('cell percent'); xticks(1:4); xtickangle(45);

plot([1 2],SOMfplc(:,PltOrd([1 2]))','-','color',[.2 .2 .2],'LineWidth',1);
plot([3 4],SOMfplc(:,PltOrd([3 4]))','-','color',[.2 .2 .2],'LineWidth',1);

title('DREADD SOM');
set(gca,'xticklabel',labels(PltOrd),'Fontsize',20,'fontname','arial','Fontweight','normal',...
    'TickDir','out','lineWidth',2,'box','off');


% STATS
fprintf('Testing %s vs. %s\n',labels{cats(1)},labels{cats(2)})
[~,pf] = fishertest(PVcontMtrx);
[~,ptp] = ttest(PVfplc(:,cats(1)),PVfplc(:,cats(2)));
[psp,~] = signrank(PVfplc(:,cats(1)),PVfplc(:,cats(2)));
[~,pta] = ttest(PVfact(:,cats(1)),PVfact(:,cats(2)));
[psa,~] = signrank(PVfact(:,cats(1)),PVfact(:,cats(2)));
fprintf('PV Fisher p = %.6f\n',pf)
fprintf('PV Placefield fraction ttest p = %.6f signrank p = %.6f\n',ptp, psp)
fprintf('PV Active fraction ttest p = %.6f signrank p = %.6f\n',pta, psa)

[~,pf] = fishertest(SOMcontMtrx);
[~,ptp] = ttest(SOMfplc(:,cats(1)),SOMfplc(:,cats(2)));
[psp,~] = signrank(SOMfplc(:,cats(1)),SOMfplc(:,cats(2)));
[~,pta] = ttest(SOMfact(:,cats(1)),SOMfact(:,cats(2)));
[psa,~] = signrank(SOMfact(:,cats(1)),SOMfact(:,cats(2)));
fprintf('SOM Fisher p = %.6f\n',pf)
fprintf('SOM Placefield fraction ttest p = %.6f signrank p = %.6f\n',ptp, psp)
fprintf('SOM Active fraction ttest p = %.6f signrank p = %.6f\n',pta, psa)

clear PVnplc PVntot PVfplc SOMnplc SOMntot SOMfplc cats PVcontMtrx SOMcontMtrx...
    PVnact PVfact SOMnact SOMfact colors PltOrd labels nSOMcells nPVcells
 
%% Axon: Get dF/F for all cells in dataset
allDF = [];
outlier_thr = 3; % 3 extreme, 1.5 moderate, Inf to disable
for ds = 1:length(AxonSomHilus)
    allDF = cat(1,allDF,AxonSomHilus{ds}.dF_F);
end

[~,foutl] = tukeyOutlierRemoval(allDF(:,1),outlier_thr);
[~,noutl] = tukeyOutlierRemoval(allDF(:,2),outlier_thr);
notOutl = find(~foutl&~noutl);
allDF = allDF(notOutl,:);

figure
distributionPlot(allDF(:,1),'histOri','left','color',[1 .7 .7],'widthDiv',[2 1],'showMM',6);
distributionPlot(allDF(:,2),'histOri','right','color',[1 .5 .5],'widthDiv',[2 2],'showMM',6);
clear ds
