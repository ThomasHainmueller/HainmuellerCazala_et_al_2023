function [resulttable, cellinfo] = metastatistics(dataset,unit,parameter,varargin)
% Create a bar plot of mean+/-SEM for a selected parameter from several
% 'tdata' type datasets collated in one cell array (metadataset) across
% the categories present in these datasets.
% units: animal, dataset, cell, hybrid
% paramters: *any field-subfield combination*, special cases:
% 'ncells' = acc to specified parameters, 'placefield_width' in m
% resulttable: individual values for all observation-units
% 
% 230806 - added 'cellinfo' to output, which gives seesion ID, dataset ID, and animal

args=struct('active',0,'SI',1.0,'PF',1.0,'region','',...
    'categories',1:length(dataset{1}.metadata.categories),'combination','or',...
    'rescategories',1:4,'plotting',true); 
    % SI/PF significance levels of spatial info/placefields from bootstrap
    % rescategories: categories for which results are returned

% Overwrite default parameters if required
for pair = reshape(varargin,2,[])
    %args.(pair{1}) = pair{2};
    if isfield(args,pair{1})
        args.(pair{1})=pair{2};
    else
        error('Input argument does not exist\n')
    end
end
resulttable=[];

% Use recursion to get cells from each animal
if strcmp(unit,'animal')
    for d=1:length(dataset)
        animalcat(d) = dataset{d}.animal;
    end
    animals = unique(animalcat);
    resultcell = {};
    % Get cellular parameters for all subdatasets acquired in this animal
    for a = length(animals):-1:1
        resultcell{a} = metastatistics(dataset(animalcat == animals(a)),...
            'cell',parameter,'active',args.active,'SI',args.SI,'PF',args.PF,...
            'region',args.region,'categories',args.categories,'combination',...
            args.combination,'rescategories',args.rescategories,'plotting',false);
        resulttable(a,args.rescategories) = nanmean(resultcell{a});
    end
end

for d=1:length(dataset)
    if strcmp(unit,'dataset')    
        for c = args.categories
            if strcmp(parameter,'ncells')
                resulttable(d,c)=length(findcells(...
                    dataset{d}, c, args.active, args.SI, args.PF, args.combination));
                %resulttable(d,c)=length(findcells(...
                %dataset{d},args.categories,args.active,args.SI,args.PF,args.combination));
            end
            
            if strcmp(parameter,'Fcells')
                resulttable(d,c)=length(findcells(...
                    dataset{d},args.categories,args.active,args.SI,args.PF,args.combination))/...
                    length(dataset{d}.cells);
            end   
            
            % nplacecells and nSI cells are redundant with the above, keep
            % for leagcy reasons for now.
            if strcmp(parameter,'nplacecells')
                resulttable(d,c)=length(find(....
                    dataset{d}.metadata.Placefield_P(:,c)<=0.05));
            end
            if strcmp(parameter,'nSIcells')
                resulttable(d,c)=length(find(....
                    dataset{d}.metadata.spatial_P(:,c)<=0.05));
            end   
        end
    elseif strcmp(unit,'cell')
        % 221210 - exclude datasets that don't have all the rescategories
        catexist = 1:length(dataset{d}.metadata.categories);
        
        % 230806 - get animal ID
        if isfield(dataset{d},'animal')
            animal = dataset{d}.animal;
        elseif isfield(dataset{d}.metadata.categories{1},'filename')
            animalstr = dataset{d}.metadata.categories{1}.filename{1};
            animalstr = strsplit(animalstr,'_');
            animal = str2num(animalstr{2});
        else
            animal = NaN;
        end
        
        if length(intersect(catexist,args.rescategories)) ~= length(args.rescategories)
            selected = [];
            
        % APPENDED none category 190911
        elseif strcmp(args.combination,'none')
            selected = 1:length(dataset{d}.cells);
            
        else
            selected=findcells(dataset{d},args.categories,args.active,args.SI,...
                args.PF,args.combination, args.region);
        end
        
        for n=1:length(selected)
            index=size(resulttable,1)+1;
            if strcmp(parameter,'Pearson_r') || strcmp(parameter,'Random_Pearson_r')
                thisres = getfield(dataset{d}.cells{selected(n)},parameter);
                resulttable(index,1:length(args.rescategories),...
                    1:length(args.rescategories)) =...
                    thisres(args.rescategories,args.rescategories);
            elseif strcmp(parameter,'placefield_width')
                % placefield width in BINS! e.g. 1=5cm with 80 bins on 4m
                for c = args.rescategories
                    resulttable(index,c) = length(find(...
                        dataset{d}.cells{selected(n)}.categories{c}.placefield));
                end
            elseif isfield(dataset{d}.cells{selected(n)},parameter)
            	result = getfield(dataset{d}.cells{selected(n)},parameter);
                resulttable(index,args.rescategories)=result(args.rescategories);
            end
            
            cellinfo(index).sessionID = selected(n);
            cellinfo(index).session = d;
            cellinfo(index).animal = animal;
        end
    elseif strcmp(unit,'hybrid') % plot cellular property averages over datasets
        thisres = []; % Empty variable from last loop.
        % APPENDED none category 190911
        if ~strcmp(args.combination,'none')
            selected=findcells(dataset{d},args.categories,args.active,args.SI,...
                args.PF,args.combination, args.region);
        else
            selected = 1:length(dataset{d}.cells);
        end
        for n=1:length(selected)
            if strcmp(parameter,'Pearson_r') || strcmp(parameter,'Random_Pearson_r')
                pearson=getfield(dataset{d}.cells{selected(n)},parameter);
                thisres(n,1:length(args.rescategories),...
                    1:length(args.rescategories)) =...
                    pearson(args.rescategories, args.rescategories);
            elseif isfield(dataset{d}.cells{selected(n)},parameter)
                param = getfield(dataset{d}.cells{selected(n)},parameter);
                thisres(n,1:length(args.rescategories))=...
                    param(1:length(args.rescategories));
            end
        end
        
        if strcmp(parameter,'Pearson_r') || strcmp(parameter,'Random_Pearson_r')
            resulttable(d,:,:) = nanmean(thisres,1);
        else
            %resulttable(d,args.categories) = nanmean(thisres,1); %changed
            %to rescategories 190821
            resulttable(d,args.rescategories) = nanmean(thisres,1);
        end
    end
end

% Remove empty columns; Pulled out of the plotting section 190809
% Found commented out on 221206, reactivated on this day.
resulttable(:,all(resulttable==0,1))=[];

% Plotting part; conditional added 190809
if args.plotting 
    try
        if strcmp(parameter,'Pearson_r') || strcmp(parameter,'Random_Pearson_r')
            figure;
            imagesc(reshape(nanmean(resulttable,1),size(resulttable,2),size(resulttable,3)));
            caxis([0,.8]);
            colormap('jet');
        else
            figure; hold on;
            bar(nanmean(resulttable,1),'b');
            errorbar(nanmean(resulttable,1),nanstd(resulttable,1)/sqrt(size(resulttable,1)),'k.');
            ylim([0,1.2*max(nanmean(resulttable,1))]);
            
            figure;
            boxplot(resulttable);
        end
    catch
        fprintf('Results not plottable\n');
    end
end
end

