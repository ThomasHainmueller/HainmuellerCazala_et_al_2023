%% Add animal ID on to all metadatasets
names = {'CA1PV','CA1SOM','CA3PV','CA3SOM','DGPV','DGSOM','DREADDpv','DREADDsom'};

for nmds = 1:length(names)
    eval(sprintf('mds = %s;',names{nmds}))
    [restable{nmds}, mds] = metadataset_sessionInfo(mds);
    eval(sprintf('%s = mds;',names{nmds})); % Overwrite old dataset var
end

clear nmds mds names

%% Get bg-corrected data from current folder
fn = dir('*_swf.mat');
load(fn.name,'data');
load('background.mat','bg');

data = correct_time(data);
data = subtract_background(data,bg);
data = split_famnov(data);

save(strrep(fn.name,'_swf.mat','_data.mat'),'data');
clear fn bg

%% Preprocess new CA1 PV datasets (210122 - individual)
fn = dir('*_swf.mat');
load(fn.name,'data');
%load('background.mat','bg');

data = correct_time(data);

for tr = 1:length(data.metadata.categories{1}.y)
    close all
    figure; hold on;
    plot(data.metadata.categories{1}.y{tr});
    %plot(diff(data.metadata.categories{1}.x{tr})*1000);
    splitpoint(tr) = input('Enter split point:\n');
    data.metadata.categories{1}.x{tr}(1:splitpoint(tr)) = .1;
    data.metadata.categories{1}.x{tr}(splitpoint(tr)+1:end) = 2;
end

data = split_famnov(data);

save(strrep(fn.name,'_swf.mat','_data.mat'),'data');
clear fn

%% Manually determine A-B transition (for new CA1 DS)
for tr = 1:length(data.metadata.categories{1}.y)
    close all
%     figure; hold on;
%     plot(data.metadata.categories{1}.y{tr});
%     %plot(diff(data.metadata.categories{1}.x{tr})*1000);
%     splitpoint(tr) = input('Enter split point:\n');
    data.metadata.categories{1}.x{tr}(1:splitpoint(tr)) = .1;
    data.metadata.categories{1}.x{tr}(splitpoint(tr)+1:end) = 2;
end

%% Homogenize new CA1 PV datasets (210122) to comply with previous data standard
for ds = 1:length(CA1PV)
    for c = 1:length(CA1PV{ds}.metadata.categories)
        if iscell(CA1PV{ds}.metadata.categories{c}.acquisition_rate)
            CA1PV{ds}.metadata.categories{c}.acquisition_rate =...
                [CA1PV{ds}.metadata.categories{c}.acquisition_rate{:}];
        end
    end
end

clear c ds

%% Convert F0 normalization to Sheffield-type
names = {'CA1PV','CA1SOM','CA3PV','CA3SOM','DGPV','DGSOM'}; %{'CA1PVnew'};

for nmds = 1:length(names)
    eval(sprintf('mds = %s;',names{nmds}))
    for nds = 1:length(mds)
        for n = 1:length(mds{nds}.cells)
            for r = 1:length(mds{nds}.cells{n}.categories{1}.dFoT)
                % Global normalization & zscoring for cat 1-2 to avoid fam/nov problems!
                oF0 = mds{nds}.cells{n}.categories{1}.F0{r}(1);
                
                thissignal = cat(2,mds{nds}.cells{n}.categories{1}.dFoT{r},...
                    mds{nds}.cells{n}.categories{2}.dFoT{r}) * oF0 + oF0;
                
                nF0 = prctile(thissignal,8);
                %Rsigma = nanstd(thissignal);
                
                thissignal = (thissignal-nF0)/nF0;
                Zsigma = nanstd(thissignal);
                Zmean = nanmean(thissignal);
                
                mds{nds}.cells{n}.categories{1}.dFoT{r} = ...
                    mds{nds}.cells{n}.categories{1}.dFoT{r} * oF0/nF0 + oF0/nF0 -1;
                mds{nds}.cells{n}.categories{1}.zscored{r} = ...
                    (mds{nds}.cells{n}.categories{1}.dFoT{r} - Zmean)/Zsigma;
                mds{nds}.cells{n}.categories{1}.F0{r} = nF0;
                %mds{nds}.cells{n}.categories{1}.baselineSD(r) = Rsigma;
                
                mds{nds}.cells{n}.categories{2}.dFoT{r} = ...
                    mds{nds}.cells{n}.categories{2}.dFoT{r} * oF0/nF0 + oF0/nF0 -1;
                mds{nds}.cells{n}.categories{2}.zscored{r} = ...
                    (mds{nds}.cells{n}.categories{2}.dFoT{r} - Zmean)/Zsigma;
                mds{nds}.cells{n}.categories{2}.F0{r} = nF0;
                %mds{nds}.cells{n}.categories{2}.baselineSD(r) = Rsigma;
            end
            
            % 'Normal' z-scoring for remaining categories
            for c = 3:length(mds{nds}.cells{n}.categories)
                for r = 1:length(mds{nds}.cells{n}.categories{c}.dFoT)
                    mds{nds}.cells{n}.categories{c}.zscored{r} = ...
                        zscore(mds{nds}.cells{n}.categories{c}.dFoT{r});
                end
            end
        end
    end
    eval(sprintf('%s = mds;',names{nmds})); % Overwrite old dataset var
end

clear nF0 oF0 sigma nds nmds n c Zsigma Zmean thissignal Rsigma names r mds

%% TBD --------- Sheffield normalization for CA1 PYR
names = {'CA1'};

for nmds = 1:length(names)
    eval(sprintf('mds = %s;',names{nmds}))
    for nds = 1:length(mds)
        for n = 1:length(mds{nds}.cells)
            for c = 1:length(mds{nds}.cells{n}.categories)
                for r = 1:length(mds{nds}.cells{n}.categories{c}.dFoT)
                    % Global normalization & zscoring for cat 1-2 to avoid fam/nov problems!
                    oF0 = mds{nds}.cells{n}.categories{c}.F0{r}(1);
                    
                    thissignal = cat(2,mds{nds}.cells{n}.categories{c}.dFoT{r}) * oF0 + oF0;
                    
                    nF0 = prctile(thissignal,8);
                    %Rsigma = nanstd(thissignal);
                    
                    thissignal = (thissignal-nF0)/nF0;
                    Zsigma = nanstd(thissignal);
                    Zmean = nanmean(thissignal);
                    
                    mds{nds}.cells{n}.categories{c}.dFoT{r} = ...
                        mds{nds}.cells{n}.categories{c}.dFoT{r} * oF0/nF0 + oF0/nF0 -1;
                    mds{nds}.cells{n}.categories{c}.zscored{r} = ...
                        (mds{nds}.cells{n}.categories{c}.dFoT{r} - Zmean)/Zsigma;
                    mds{nds}.cells{n}.categories{c}.F0{r} = nF0;
                    %mds{nds}.cells{n}.categories{1}.baselineSD(r) = Rsigma;
                end
            end
        end
    end
    eval(sprintf('%s = mds;',names{nmds})); % Overwrite old dataset var
end

clear nF0 oF0 sigma nds nmds n c Zsigma Zmean thissignal Rsigma names r mds

%% Correct moving periods
names = {'CA1PV','CA1SOM','CA3PV','CA3SOM','DGPV','DGSOM'};%{'CA1PVnew'};%
minspeed = 2;

for nmds = 1:length(names)
    eval(sprintf('mds = %s;',names{nmds}))
    for nds = 1:length(mds)
        for c = 1:length(mds{nds}.metadata.categories)
            for r = 1:length(mds{nds}.metadata.categories{c}.y)
                fr = mds{nds}.metadata.categories{c}.acquisition_rate(r);
                thisy = mds{nds}.metadata.categories{c}.y{r}*200;
                thisv = movmedian([0 abs(diff(movmedian(thisy,round(fr/3))))]*fr,fr);
                mds{nds}.metadata.categories{c}.moving{r} = thisv > minspeed;
            end
        end
    end
    eval(sprintf('%s = mds;',names{nmds}))
end
clear thisy thisv mds nds nmds c r minspeed names

%% Update datasets with transientrate, spatial info, etc.
names = {'CA1PV','CA1SOM','CA3PV','CA3SOM','DGPV','DGSOM'};%{'CA1'};

for nmds = 1:length(names)
    eval(sprintf('mds = %s;',names{nmds}))
    for nds = 1:length(mds)
        % fprintf('Recalculating %s - dataset #%i\n',names{nmds},nds)
        for n = 1:length(mds{nds}.cells)
            for c = 1:length(mds{nds}.cells{n}.categories)
                for r = 1:length(mds{nds}.metadata.categories{c}.y)
                    mds{nds}.cells{n}.categories{c}.transientmask{r}=...
                        true(size(mds{nds}.metadata.categories{c}.y{r}));
                end
                signals = cat(2, mds{nds}.cells{n}.categories{c}.dFoT{:});       
                mds{nds}.cells{n}.AUCrate(c) = sum(signals) / length(signals) *...
                    mean(mds{nds}.metadata.categories{c}.acquisition_rate);
                mds{nds}.cells{n}.transientrate(c) = 0;
                mds{nds}.cells{n}.Placefield_P(c) = 1;
            end
        end
        
        mds{nds} = update_dataset(mds{nds});
        mds{nds} = recalculate_dFoY(mds{nds},1,0,0.1:0.025:1.9,0);
        % mds{nds} = spatial_info_bootstrap(mds{nds},1000,50,0.1:0.025:1.9);
        % NOTE: SI and bootstrap are calculated later (circular shuffle)

%         for n=1:length(mds{nds}.cells)
%             mds{nds}.cells{n}.SIperAUC = ...
%                 mds{nds}.cells{n}.spatialinfo./mds{nds}.cells{n}.AUCrate;
%         end
        eval(sprintf('%s = mds;',names{nmds}))
    end
end

clear c n r nds nmds signals names fr mds

%% Running-speed profiles
names = {'CA1PV','CA1SOM','CA3PV','CA3SOM','DGPV','DGSOM'};%{'CA1'};%

for nmds = 1:length(names)
    eval(sprintf('mds = %s;',names{nmds}))
    for nds = 1:length(mds)
        mds{nds} = recalculate_dFodY(mds{nds}); % Add dF/dZ over speed profiles
        mds{nds} = recalculate_dFoY(mds{nds},true,false,.1:.025:1.9,true); % Add dZoY profiles!; Param data,moving,transients,bins,zscored
    end
    eval(sprintf('%s = mds;',names{nmds}))
end
clear mds nds nmds c r names

%% Speed modulation
names = {'CA1PV','CA1SOM','CA3PV','CA3SOM','DGPV','DGSOM'};
xv = (0:15)';

for nmds = 1:length(names)
    eval(sprintf('mds = %s;',names{nmds}))

    for nds = 1:length(mds)
        for n = 1:length(mds{nds}.cells)
            for c = 1:length(mds{nds}.metadata.categories)
                ftr = nanmean(cat(2,...
                    mds{nds}.cells{n}.categories{c}.dFodY{:}),2);
                ztr = nanmean(cat(2,...
                    mds{nds}.cells{n}.categories{c}.dZodY{:}),2);
                ftr = ftr(1:16); ztr = ztr(1:16);
                
                ftr = movmean(ftr,3,'omitnan');
                ztr = movmean(ztr,3,'omitnan');
                fnan = isnan(ftr); znan = isnan(ztr);
                %ftr(fnan) = 0; ztr(znan) = 0;
                
                Pfit = polyfit(xv(~fnan),ftr(~fnan),1);
                mds{nds}.cells{n}.speedmod(c) = Pfit(1);
                mds{nds}.cells{n}.speedYincept(c) = Pfit(2);
                
                % Correlation
                [R,p] = corrcoef(xv(~fnan),ftr(~fnan));
                mds{nds}.cells{n}.speed_R(c) = R(1,2);
                mds{nds}.cells{n}.speed_P(c) = p(1,2);
            end
        end
    end
    
    % WRITE BACK TO ORIGINAL DATA (.speedmod and .speedYintercep,
    % .speed_R,. speed_P are added to cell)
    eval(sprintf('%s = mds;',names{nmds}))
end

clear mds nds nmds n c r names DF DZ DFc gn fi fnan ftr P xv yv znan ztr p R Pfit

%% Mobile - immobile ratio
names = {'CA1PV','CA1SOM','CA3PV','CA3SOM','DGPV','DGSOM'};

for nmds = 1:length(names)
    eval(sprintf('mds = %s;',names{nmds}))
    for nds = 1:length(mds)
        for c = 1:length(mds{nds}.metadata.categories)
            allmoving = cat(2,mds{nds}.metadata.categories{c}.moving{:});
            for n = 1:length(mds{nds}.cells)
                allsignals = cat(2,mds{nds}.cells{n}.categories{c}.dFoT{:});
                
                Amov = nanmean(allsignals(allmoving));
                Aimm = nanmean(allsignals(~allmoving));
                mds{nds}.cells{n}.MovImmRatio(c) = Amov/Aimm;
                mds{nds}.cells{n}.Amov(c) = Amov;
                mds{nds}.cells{n}.Aimm(c) = Aimm;
            end
        end
    end
    
    % SAVE back to dataset in workspace, omit if you just want stats
    %eval(sprintf('%s = mds;',names{nmds}))
    
    allmv = metastatistics(mds,'cell','Amov','plotting',false,'categories',1,'rescategories',1);
    allim = metastatistics(mds,'cell','Aimm','plotting',false,'categories',1,'rescategories',1);
    p(nmds) = signrank(allmv,allim);
    values{nmds} = allmv./allim;
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

clear allmoving allsignals Amov Aimm mds nmds nds c n names allv values allim allmv p

%% Spatial coherence, etc.
names = {'CA1PV','CA1SOM','CA3PV','CA3SOM','DGPV','DGSOM'};

for nmds = 1:length(names)
    eval(sprintf('mds = %s;',names{nmds}))
    for nds = 1:length(mds)
        mds{nds} = spatial_coherence(mds{nds});
        mds{nds} = reviewer2_corr(mds{nds});
    end
    eval(sprintf('%s = mds;',names{nmds}))
end
clear names nmds nds mds

%% Add Olypher positional information score
names = {'CA1PV','CA1SOM','CA3PV','CA3SOM','DGPV','DGSOM','CA1'};
for nmds = 1:length(names)
    eval(sprintf('mds = %s;',names{nmds}))
    for nds = 1:length(mds)
        fprintf('Processing %s dataset# %d\n',names{nmds},nds);
        mds{nds} = PositionalInfo(mds{nds});
    end
    eval(sprintf('%s = mds;',names{nmds}))
end

clear mds nmds nds names

%% Calculate pos and neg place fields (reviewer 3)
names = {'CA1PV','CA1SOM','CA3PV','CA3SOM','DGPV','DGSOM','CA1'};
for nmds = 1:length(names)
    eval(sprintf('mds = %s;',names{nmds}))
    for nds = 1:length(mds)
        fprintf('Processing %s dataset# %d\n',names{nmds},nds);
        mds{nds} = PlacefieldsPosNeg(mds{nds});
    end
    eval(sprintf('%s = mds;',names{nmds}))
end

clear mds nmds nds names

%% Calculate activity-difference-scores and between-context-correlations (remapping)
names = {'CA1PV','CA1SOM','CA3PV','CA3SOM','DGPV','DGSOM'};

for nmds = 1:length(names)
    eval(sprintf('mds = %s;',names{nmds}))
    for nds = 1:length(mds)
        % Correlations of spatial maps between categories
        mds{nds} = spatial_corr(mds{nds},1,0,1,1);
        
        % Activity difference scores
        for n = 1:length(mds{nds}.cells)
            mds{nds}.cells{n}.AUCdiffScore = ...
                abs((mds{nds}.cells{n}.AUCrate(1) - mds{nds}.cells{n}.AUCrate(2))/...
                (mds{nds}.cells{n}.AUCrate(1) + mds{nds}.cells{n}.AUCrate(2)));
            
            % F-N correlation
            mds{nds}.cells{n}.FNcorr = mds{nds}.cells{n}.Pearson_r(1,2);
        end
    end
    eval(sprintf('%s = mds;',names{nmds}))
end
clear names nmds nds mds

%% Spatial parameters for randomly shifted traces
% Moving/Immobile ratio, speed-modulation slope, SI, SI/AUC, trial-trial
% variability, 1st to 2nd half stability
names = {'CA1PV','CA1SOM','CA3PV','CA3SOM','DGPV','DGSOM'};

for nmds = 1:length(names)
    eval(sprintf('mds = %s;',names{nmds}))
    for nds = 1:length(mds)
        fprintf('Processing dataset# %d\n',nds);
        mds{nds} = spatial_param_rnd(mds{nds},'category',1,'nshuffles',1000);
    end
    eval(sprintf('%s = mds;',names{nmds}))
end
clear names nmds nds mds

%% Make random parameters useful for metastatistics
names = {'CA1PV','CA1SOM','CA3PV','CA3SOM','DGPV','DGSOM'};
% names = {'DGSOM'};

for nmds = 1:length(names)
    eval(sprintf('mds = %s;',names{nmds}))
    for nds = 1:length(mds)
        for n = 1:length(mds{nds}.cells)
            mds{nds}.cells{n}.spatialinfo_r(1) = ...
                mds{nds}.cells{n}.spatialinfo_rand{1}(1);
            mds{nds}.cells{n}.SIperAUC(1) = ...
                mds{nds}.cells{n}.spatialinfo(1) / mds{nds}.cells{n}.AUCrate(1);
            mds{nds}.cells{n}.SIperAUC_r(1) = ...
                mds{nds}.cells{n}.spatialinfo_r / mds{nds}.cells{n}.AUCrate(1);
            mds{nds}.cells{n}.spatial_coherence_r(1) = ...
                mds{nds}.cells{n}.spatial_coherence_rand{1}(1);
            mds{nds}.cells{n}.sessionstab_r(1) = ...
                mds{nds}.cells{n}.sessionstab_rand{1}(1);
            mds{nds}.cells{n}.tvl_r(1) = ...
                mds{nds}.cells{n}.TVLr{1}(1);
        end
    end
    eval(sprintf('%s = mds;',names{nmds}))
end
clear names nmds nds mds

%% Get cell numbers and animal per dataset
names = {'CA1PV','CA1SOM','CA3PV','CA3SOM','DGPV','DGSOM'};
%f1 = figure('Position',[200,200,600,600]);
[PVanimals, SOManimals] = deal({});

for nmds = 1:length(names)
    eval(sprintf('mds = %s;',names{nmds}))
    for nds = 1:length(mds)
        ncells(nmds,nds) = length(mds{nds}.cells);
        %fprintf('%s dataset #%i name = %s\n',names{nmds},nds,mds{nds}.metadata.categories{1}.filename{1});
        nm = mds{nds}.metadata.categories{1}.filename{1};
        nm = strsplit(nm,'_');
        animal{nmds,nds} = nm{2};
        
        if mod(nmds,2)
            PVanimals{end+1} = nm{2};
        else
            SOManimals{end+1} = nm{2};
        end
        
        if ~isfield(mds{nds},'animal')
            mds{nds}.animal = animal;
        end
    end
   
    %eval(sprintf('%s = mds;',names{nmds}))
end

fprintf('Mean cell number: %.2f, range %i to %i\n',mean(ncells,'all'),...
    min(ncells(ncells>0)),max(ncells,[],'all'));

fprintf('Total animal count: %i PV-mice and %i SOM mice\n',length(unique(PVanimals)),...
    length(unique(SOManimals)));

clear names nmds nds nm

%% Get session and cell number organized by animal
dsname = 'DREADDctr';
region = 'DG';

eval(sprintf('mds = %s;',dsname));
tab = metadataset_animalInfo(mds);

for r = 1:length(tab)
    fprintf('%i\t%s (%i)\t%s (%i)\n',...
        tab(r).animal,region,tab(r).nsessions,region,tab(r).ncells)
end

clear mds dsname region tab

%% DREADD IN series: Update AUC rate following df_f_meta_soma_ac
dsnames = {'DREADDinPV','DREADDinSOM'};
for nmds = 1:length(dsnames)
    eval(sprintf('mds = %s;',dsnames{nmds}));
    for ds = 1:length(mds)
        for c = 1:length(mds{ds}.metadata.categories)
            moving = cat(2,mds{ds}.metadata.categories{c}.moving{:});
            
            for n = 1:length(mds{ds}.cells)
                trace = cat(2,mds{ds}.cells{n}.categories{c}.dFoT_fneu{:});
                trace(trace<0) = 0;
                mds{ds}.cells{n}.AUCrate(c) = nanmean(trace(moving));
            end
        end
    end
    eval(sprintf('%s = mds;',dsnames{nmds}));
end

clear mds nmds ds n c moving trace dsnames

%% Kruskal-Wallis with Dunn's
label = {'CA1PV','CA1SOM','CA3PV','CA3SOM','DGPV','DGSOM'}; %'CA1PYR',
valtbl = catuneven(values,NaN);
[p,tbl,stats] = kruskalwallis(valtbl);
[compc,m] = multcompare(stats,'CType','dunn-sidak');

% Output labelled group comparison table!
fprintf('\nNumber of observations (N)\n')
for n = 1:length(label)
    fprintf('%s: %i\n',label{n},stats.n(n))
end

fprintf('\nKruskal-Wallis\n')
fprintf('%s\t%s\t%s\t%s\t%s\t%s\n',tbl{1,1},tbl{1,2},tbl{1,3},tbl{1,4},tbl{1,5},tbl{1,6})
for n = 2:size(tbl,1)
    fprintf('%s\t%.2e\t%i\t%.2e\t%.2f\t%.2e\n',tbl{n,1},tbl{n,2},tbl{n,3},tbl{n,4},tbl{n,5},tbl{n,6})
end

fprintf('\nDunn\''s post-hoc test\n')
fprintf('%s\t%s\t%s\n','Group1','Group2','p')
for n = 1:size(compc,1)
    fprintf('%s\t%s\t%.2e\n',label{compc(n,1)},label{compc(n,2)},compc(n,6))
end

clear label n

%% OBSOLETE ------- Recalculate factors for Multifactorial analysis
% Include:
% (3) tuning-vector length, (4) session stability (6) spatial coherence
names = {'CA1PV','CA1SOM','CA3PV','CA3SOM','DGPV','DGSOM'};

for nmds = 1:length(names)
    eval(sprintf('mds = %s;',names{nmds}))
    
    
    
    for p = 1:length(parameters)
        results{nmds} = cat(2,results{nmds},metastatistics(mds,'cell',...
            parameters{p},'categories',1,'rescategories',1,'plotting',false));
    end
    %eval(sprintf('%s = mds;',names{nmds}))
end

clear names nmds mds

%% Multifactorial analysis and parameter correlation for interneurons
% Run 'IN_plots.m', section 'FamNov Comparison' before this!
%
% Include:
% (1) Mov/Imm ratio, (2) speed modulation slope
% (3) tuning-vector length, (4) session stability (5) trial-trial-corr 
% (6) spatial coherence (7) spatial info (8) F/N ratio
% Not included: 
% - FN correlation (poor spatial tuning to begin with)
% - Normalized SI (no other activity measures included)
% - ZdiffScore (tried, conceptually better than FNratio but results are
% mostly same)

names = {'CA1PV','CA1SOM','CA3PV','CA3SOM','DGPV','DGSOM'};
colors = [0 .5 0; .5 1 .5; 0 0 1; .5 .5 1; 1 0 0; 1 .5 .5];
parameters = {'MovImmRatio','speedmod','spatialinfo','tvl','spatial_coherence',...
    'sessionstab','trialvar','FN_ratio'};
paramnames = {'MvIm','VelM','SI','Tvl','SC','Sstab','Ttvar','FNr'};

f1 = figure;
f2 = figure;
resIdcs = [];

for nmds = 1:length(names)
    eval(sprintf('mds = %s;',names{nmds}))
    results{nmds} = [];
    
    for p = 1:length(parameters)
        results{nmds} = cat(2,results{nmds},metastatistics(mds,'cell',...
            parameters{p},'categories',1,'rescategories',1,'plotting',false));
    end
    
    % Remove extreme values >3 interquartiles
    results{nmds} = tukeyOutlierRemoval(results{nmds},3);
    resIdcs(end+1:end+length(results{nmds})) = nmds;
    
    [R(:,:,nmds),pval(:,:,nmds)] = corrcoef(results{nmds});
    
    % PCA (standardized variables, equivalent to correlation matrix)
    [coeff(:,:,nmds),~,~,~,explained(:,nmds)] = pca(zscore(results{nmds}));
    
    % Plot correlations
    figure(f1);
    subplot(3,2,nmds)
    imagesc(R(:,:,nmds)); colormap(bwr); caxis([-1 1]);
    xticks(1:length(parameters)); xticklabels(paramnames); xtickangle(45);
    yticks(1:length(parameters)); yticklabels(paramnames); ytickangle(45);
    title(names{nmds});
    
    % Plot Principal components and factor loadings
    figure(f2);
    subplot(3,2,1); hold on;
    plot(1:length(parameters),explained(:,nmds),'o-','color',colors(nmds,:));
    ylabel('Explained variance (%)');
    xlabel('Principal component #');
    
    for pc = 1:5
        subplot(3,2,pc+1); hold on;
        title(sprintf('PC #%i',pc));
        plot(1:length(parameters),coeff(:,pc,nmds),'o-','color',colors(nmds,:)); % Signed
        %plot(1:length(parameters),abs(coeff(:,pc,nmds)),'o-','color',colors{nmds});% Absolute
        ylabel('Factor loading')
        xticks(1:length(parameters));
        xticklabels(paramnames);
        xtickangle(45);
    end
end

% Plot global averages
allres = cat(1,results{:});
[Rall,pall] = corrcoef(allres);
[CoeffAll,ScoreAll,~,~,explAll] = pca(zscore(allres));

figure(f2)
subplot(3,2,1);
plot(1:length(parameters),explAll,'o-','color','k','MarkerFaceColor','k');

for pc = 1:5
    subplot(3,2,pc+1); hold on;
    plot(1:length(parameters),CoeffAll(:,pc),'o-','color','k','MarkerFaceColor','k');
    %plot(1:length(parameters),abs(CoeffAll(:,pc)),'o-','color','k','MarkerFaceColor','k');
end

% Plot correlation matrix and p-values across all interneurons
f3 = figure;
subplot(2,1,1);
imagesc(Rall(:,:)); colormap(bwr); caxis([-1 1]);
title('All Interneurons');
xticks(1:length(parameters)); xticklabels(paramnames); xtickangle(45);
yticks(1:length(parameters)); yticklabels(paramnames); ytickangle(45);

subplot(2,1,2)
imagesc(pall(:,:)); caxis([0 0.05]);
xticks(1:length(parameters)); xticklabels(paramnames); xtickangle(45);
yticks(1:length(parameters)); yticklabels(paramnames); ytickangle(45);

% Scatterplot of IN parameters in cardinal factors and PC space
f4 = figure;

subplot(2,1,1);
scatter3(ScoreAll(:,1),ScoreAll(:,2),ScoreAll(:,3),[],colors(resIdcs,:),'.');
title('PC scores');
xlabel('PC#1'); ylabel('PC#2'); zlabel('PC#3');
    
subplot(2,1,2);
scatter3(allres(:,2),allres(:,8),allres(:,4),[],colors(resIdcs,:),'.');
title('Cardinal factors');
xlabel('VelM'); ylabel('FNr'); zlabel('Tvl');
set(gca, 'Yscale', 'log');

% % Factor loadings (pictorial);
% figure(f2)
% for pc = 1:5
%     subplot(2,3,pc+1); %hold on;
%     imagesc(1:length(parameters),1:length(names),...
%         abs(squeeze(coeff(:,pc,:)))');
%     title(sprintf('PC #%i factor loadings',pc));
%     colormap('hot'); caxis([-1 1]);
%     yticks(1:length(names));
%     yticklabels(names);
%     ytickangle(45);
%     xticks(1:length(parameters));
%     xticklabels(paramnames);
%     xtickangle(45);
% end

clear names parameters nmds mds p pc paramnames Rall pall CoeffAll ...
    ScoreAll

%% DEPRECATED ---- ORIGINAL SPEED PROFILES COMBINED WITH PLOTTING
names = {'CA1PV','CA1SOM','CA3PV','CA3SOM','DGPV','DGSOM'};
f1 = figure; f2 = figure; f3 = figure;
colors = {[0 1 1],[1 1 0],[.5 .5 .5]};
%plo = [1 4 2 5 3 6];
plo = 1:6;
xv = (0:15)';

% Make color bar figure
colbar(1,:,:) = ones(256,3)-repmat(0:1/255:1,3,1)'.*colors{1};
colbar(2,:,:) = ones(256,3)-repmat(0:1/255:1,3,1)'.*colors{2};
colbar(3,:,:) = ones(256,3)-repmat(0:1/255:1,3,1)'.*colors{3};
figure;
figure; image(colbar);
box off

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
                %ftr(fnan) = 0; ztr(znan) = 0;
                
                DF(:,gn,c) = ftr;
                DZ(:,gn,c) = ztr;
                
                Pfit = polyfit(xv(~fnan),ftr(~fnan),1);
                mds{nds}.cells{n}.speedmod(c) = Pfit(1);
                mds{nds}.cells{n}.speedYincept(c) = Pfit(2);
                
                % Correlation
                [R,p] = corrcoef(xv(~fnan),ftr(~fnan));
                mds{nds}.cells{n}.speed_R(c) = R(1,2);
                mds{nds}.cells{n}.speed_P(c) = p(1,2);
                
                % Make colored plot
                ztr = ztr-min(ztr);
                ztr = ztr/max(ztr);
                
                if p(1,2) < .05 && R(1,2) > 0
                    DFc(gn,:,:,c) = ones(size(ztr,1),3)-ztr*colors{1};
                elseif p(1,2) < .05 && R(1,2) < 0
                    DFc(gn,:,:,c) = ones(size(ztr,1),3)-ztr*colors{2};
                else
                    DFc(gn,:,:,c) = ones(size(ztr,1),3)-ztr*colors{3};
                end
                
                %order{nmds}(gn,c) = R(1,2);
                order{nmds}(gn,c) = Pfit(1);
            end
        end
    end
    
    for c = 1:length(mds{nds}.metadata.categories)
%         [~,maxbin] = max(DF(1:15,:,c),[],1);
%         [~,order{nmds}(:,c)] = sort(maxbin);
        [~,order{nmds}(:,c)] = sort(order{nmds}(:,c));
    end
    
    DFplot{nmds} = DF;
    DZplot{nmds} = DZ;
    DFcplot{nmds} = DFc;
    
    yv = (1:size(DFplot{nmds},2))';
    
    figure(f1);
    subplot(3,2,plo(nmds));
    imagesc(xv,yv,DZplot{nmds}(:,order{nmds}(:,1),1)'); 
    caxis([-.5 1.5]); colormap(hot); xlim([0 15]);
    xlabel('Speed (cm/s)'); ylabel('Cell #');
    title(names{nmds});
    colorbar;
    
    figure(f2);
    subplot(3,2,plo(nmds));
    imagesc(xv,yv,DFplot{nmds}(:,order{nmds}(:,1),1)'); 
    caxis([0 .7]); xlim([0 15]);
    xlabel('Speed (cm/s)'); ylabel('Cell #');
    title(names{nmds});
    colorbar
    
    figure(f3);
    subplot(3,2,plo(nmds));
    image(xv,yv,DFcplot{nmds}(order{nmds}(:,1),:,:,1)); 
    %caxis([0 .7]); xlim([0 15]);
    xlabel('Speed (cm/s)'); ylabel('Cell #');
    title(names{nmds});
    %colorbar
    
    % WRITE BACK TO ORIGINAL DATA (.speedmod and .speedYintercep are added
    % to cell).
    eval(sprintf('%s = mds;',names{nmds}))
    
    clear DF DZ DFc
end
clear mds nds nmds n c r names DF DZ DFc gn fi fnan ftr P xv yv znan ztr ...
    colors DFcplot DFplot DZplot f1 f2 f3 order p plo R
