function dFoYplot = plot_speedmaps_meta(dataset,sortcategory,varargin)
% Plots lineartrack acquired activity grouped to category for multiple
% datasets. Takes dFoY traces as given in the dataset - should allready be
% recalculated only for the significant transients. Parameters (activity,
% SI or PF p-values) are checked for the given categories and cells
% meeting the criteria will be plotted in an 'and' or 'or' manner as
% specified. Generate appropriate input by placing several tdata
% (derived from get_lineartracksX => convert_data) into a cell array.

if nargin < 2
    sortcategory = 1;
end

args=struct('active',0,'SI',1.0,'PF',1.0,'categories',sortcategory,...
    'combination','and','signal','zscored','plotmean',true,'DispCats',0); 
    % SI/PF significance levels of spatial info/placefields from bootstrap

% Overwrite default parameters if required
for pair = reshape(varargin,2,[])
    %args.(pair{1}) = pair{2};
    if isfield(args,pair{1})
        args.(pair{1})=pair{2};
    else
        error('Input argument does not exist\n')
    end
end

dFoYplot = [];

for d=1:numel(dataset)
    ds=dataset{d};
    % Initialize logical array for selected cells
    if strcmp(args.combination,'and')
        selection=true(1,length(ds.cells));
    elseif strcmp(args.combination,'or')
        selection=false(1,length(ds.cells));
    elseif strcmp(args.combination,'none')
        selection=true(1,length(ds.cells));
    else
        error('Combination parameter must be ''and'', ''or'' or ''none''');
    end
    % Check if conditions are met in all/any of selected categories
    for c = args.categories
        if strcmp(args.combination,'and')
            [~,thissel]=findcells(ds,c,args.active,args.SI,args.PF);
            selection=selection&thissel;
        elseif strcmp(args.combination,'or')
            [~,thissel]=findcells(ds,c,args.active,args.SI,args.PF);
            selection=selection|thissel;
        end
    end
    selection=find(selection); % Convert to indices
    currentn=size(dFoYplot,2);
    
    % Actual collection of the signals
    for c=length(ds.metadata.categories):-1:1
        for n=length(selection):-1:1
            if strcmp(args.signal,'zscored') % CAVE: For PVcells only, no transients implemented yet
                thisdFodY=[];
                for run = 1:length(ds.cells{selection(n)}.categories{c}.zscored)
                    if isfield(ds.cells{selection(n)}.categories{c},'dZodY')
                        rawdFodY = ds.cells{selection(n)}.categories{c}.dZodY{run};
                    else
                        % DO NOT USE - insensitive to acquisition rate
                        rawdFodY=dFoSpeed(ds.cells{selection(n)}.categories{c}.zscored{run},...
                            ds.metadata.categories{c}.y{run}, 0:1.25E-3:4E-2); % Compute zscore over dY
                    end
                    thisdFodY=cat(2,thisdFodY,rawdFodY);
                end
            end
            
            if strcmp(args.signal,'raw')
                thisdFodY=[];
                for run = 1:length(ds.cells{selection(n)}.categories{c}.dFoY) 
                    if isfield(ds.cells{selection(n)}.categories{c},'dFodY')
                        rawdFodY = ds.cells{selection(n)}.categories{c}.dFodY{run}; % Use pre computed dFodY
                    else
                        % DO NOT USE - insensitive to acquisition rate
                        rawdFodY=dFoSpeed(ds.cells{selection(n)}.categories{c}.dFoY{run},...
                            ds.metadata.categories{c}.y{run}, 0:1.25E-3:4E-2); % Compute zscore over dY
                    end
                    thisdFodY=cat(2,thisdFodY,rawdFodY);
                end
            end            
            dFoYplot(:,currentn+n,c)=nanmean(thisdFodY,2);
        end
    end
end

% Find order of placefileds in selected category
[~, sortvector] = max(dFoYplot(1:15,:,sortcategory)); %Sorting matched to display value!
%[~, sortvector] = max(dFoYplot(:,:,sortcategory)); 
[~, sortindices] = sort(sortvector);

sumfig = figure('Position',[500,250,500,400]);

if args.DispCats == 0
    args.DispCats = 1:length(dataset{1}.metadata.categories);
end

%for c = 1:length(ds.metadata.categories)
for c = args.DispCats
    FigHandle = figure('Position',[(c-1)*450,50,500,900]);
    if args.plotmean
        subplot(2,1,1)
            imagesc(dFoYplot(:,sortindices,c)');
            xlim([1,15])
            title(['Category ' num2str(c)])
            if strcmp(args.signal,'zscored')
                caxis([0.0 1]);
                colormap(jet);
            else
                caxis([-.1 0.5]);
                colormap(jet);
                %caxis([0.0 0.35]);
                %colormap(flipud(gray));
            end
         subplot(2,1,2)
             % Mean +/- std for dF over y coord for all runs together
             plot(nanmean(dFoYplot(:,:,c),2)','k')
             hold on
             plot(nanmean(dFoYplot(:,:,c),2)'...
                 -(nanstd(dFoYplot(:,:,c)')),'r');
                %-(nanstd(dFoYplot(:,:,c)')/size(dFoYplot,2)),'r');
             hold on
             plot(nanmean(dFoYplot(:,:,c),2)'...
                 +(nanstd(dFoYplot(:,:,c)')),'r');
                %+(nanstd(dFoYplot(:,:,c)')/size(dFoYplot,2)),'r');
             ylim([-.5,2.0])
             xlim([0,16])
         figure(sumfig)
             % Mean +/- std for dF over y coord for all runs together
             colorval=c/length(args.DispCats); % Color gradient: 1st = blue, last = red!
             plot(nanmean(dFoYplot(:,:,c),2)','Color',[colorval,0,1-colorval])
             hold on
             plot(nanmean(dFoYplot(:,:,c),2)'...
                 -(nanstd(dFoYplot(:,:,c)')/sqrt(size(dFoYplot,2))),...
                 'Color',[colorval,0,1-colorval]);
             hold on
             plot(nanmean(dFoYplot(:,:,c),2)'...
                 +(nanstd(dFoYplot(:,:,c)')/sqrt(size(dFoYplot,2))),...
                 'Color',[colorval,0,1-colorval]);
             %ylim([-.1,.5])
             xlim([1,16])            
    else
    	imagesc(dFoYplot(:,sortindices,c)');
        xlim([1,30])
        title(['Category ' num2str(c)])
        if strcmp(args.signal,'zscored')
            caxis([0.0 1]);
            colormap(jet);
        else
            caxis([-.1 0.2]);
            colormap(jet);
            %colormap(flipud(gray));
        end
    end
end
end