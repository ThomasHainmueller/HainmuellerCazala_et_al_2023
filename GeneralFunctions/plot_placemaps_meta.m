function [dFoYplot, sortindices] = plot_placemaps_meta(dataset,sortcategory,varargin)
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
    'combination','and','signal','dFoY','region','','bins',[0.1:0.025:2.1],...
    'runmean',1,'clims',[],'xscalefactor',2,'figscale',[1 1]); 
    % SI/PF significance levels of spatial info/placefields from bootstrap
    % Combination: 'and', 'or', 'none' (to avoid condition checking)
    % Signal: 'dFoY (default)' or 'zscored'
    % 'xscalefactor' = ratio between x-unit and distance, default 2.0

% Overwrite default parameters if required
for pair = reshape(varargin,2,[])
    %args.(pair{1}) = pair{2};
    if isfield(args,pair{1})
        args.(pair{1})=pair{2};
    else
        error('Input argument does not exist\n')
    end
end

filter = fspecial('gaussian',[9,1],1.0);
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
    % Check if conditions are met in all/any of selected categories =>
    % TODO: May be simplified using full 'findcells' functionality!
    for c = args.categories
        if strcmp(args.combination,'and')
            [~,thissel]=findcells(ds,c,args.active,args.SI,args.PF,'and',args.region);
            selection=selection&thissel;
        elseif strcmp(args.combination,'or')
            [~,thissel]=findcells(ds,c,args.active,args.SI,args.PF,'or',args.region);
            selection=selection|thissel;
        end
    end
    selection=find(selection); % Convert to indices
    currentn=size(dFoYplot,2);
    
    for c=length(ds.metadata.categories):-1:1
        for n=length(selection):-1:1
            if strcmp(args.signal,'dFoY')
                thisdFoY=deal(ds.cells{selection(n)}.categories{c}.dFoY);
                thisdFoY=cat(2,thisdFoY{:});
            elseif strcmp(args.signal,'zscored') % CAVE: For PVcells only, no transients implemented yet
                thisdFoY=[];
                for run = 1:length(ds.cells{selection(n)}.categories{c}.zscored) 
                    rawdFoY = SBdiscretize(ds.cells{selection(n)}.categories{c}.zscored{run},...
                        ds.metadata.categories{c}.y{run},args.bins);
                    thisdFoY=cat(2,thisdFoY,rawdFoY);
                end
            end
            % Apply Gaussian filter for each trace
            for tr=1:size(thisdFoY,2)
                thisline=thisdFoY(:,tr);
                nanvals = isnan(thisline);
                thisline(nanvals)=0;
                thisline = conv(thisline,filter,'same');
                thisline(nanvals)=NaN;
                thisdFoY(:,tr)=thisline;
            end
            dFoYplot(:,currentn+n,c)=nanmean(thisdFoY,2);
        end
    end
end

% Smooth output if desired
dFoYplot = movmean(dFoYplot,args.runmean,1,'omitnan');

% Find order of placefileds in selected category
[~, sortvector] = max(dFoYplot(:,:,sortcategory));
[~, sortindices] = sort(sortvector);

for c = 1:length(ds.metadata.categories)
     FigHandle = figure('Position',[(c-1)*320,300,300*args.figscale(1),500*args.figscale(2)]);
     
     % Modified 221209 for nicer axis labels, etc.
     thisplot = dFoYplot(:,sortindices,c)';
     
     if ~isempty(args.xscalefactor)
         % Make sure values start at zero
         xv = (args.bins - args.bins(1)) * args.xscalefactor;
     else
         xv = 1:size(thisplot,2);
     end
     yv = 1:size(thisplot,1);
     
     imagesc(xv,yv,thisplot);
     title(['Category ' num2str(c)])
     
     xlabel('Track distance (m)');
     yticklabels('');
     unitstr = sprintf('N = %i',yv(end));
     text(0.3,0.95*yv(end),unitstr,'color',[1 1 1],'FontSize',10)
     
     if ~isempty(args.clims) % Introduced 221209
         caxis(args.clims);
         colormap(jet);
     elseif strcmp(args.signal,'zscored')
         %caxis([0.0 1.5]);
         caxis([0,1.3]);
         %colormap(flipud(gray));
         colormap('hot');
     else
        caxis([0.0 0.35]);
        %colormap(flipud(gray));
        colormap(jet);
     end
end
end
 
% for c = length(data.metadata.categories):-1:1
%     cats{c}.alltracks = [];
%     newcellno = 1;
%     % Loop over cells and collect activities over time and Y
%     for cell = find(Lindices)
%         thisdFoT = [];
%         thisdFoY = [];
%         thismask = [];
%         for run = 1:length(data.cells{cell}.categories{c}.dFoY) 
%             if transients
%                 thismask = data.cells{cell}.categories{c}.transientmask{run}';
%                 thisdFoT = cat(1,thisdFoT,data.cells{cell}.categories{c}.dFoT{run}'.*thismask);
%                 rawdFoY = data.cells{cell}.categories{c}.dFoY{run};
%             else
%                 thisdFoT = cat(1,thisdFoT,data.cells{cell}.categories{c}.dFoT{run}');
%                 rawdFoY = discretize(data.cells{cell}.categories{c}.dFoT{run},...
%                     data.metadata.categories{c}.y{run},0.1:0.05:2.1);
%             end
%             % Use Gaussian filter to smooth the space plot.
%             nanvals = isnan(rawdFoY);
%             rawdFoY(nanvals)=0;
%             rawdFoY = conv(rawdFoY,filter,'same');
%             rawdFoY(nanvals)=NaN;
%             thisdFoY = cat(2,thisdFoY,rawdFoY);
%         end
%         % Append dFoT activity from this cell for all trials.
%         cats{c}.alltracks = cat(2,cats{c}.alltracks,thisdFoT);
%         
%         % Average dFoY activity over all trials
%         cats{c}.spaceplot(:,newcellno) = nanmean(thisdFoY,2);
%         newcellno=newcellno+1;
%         clear thisdFoT thisdFoY
%     end
% end
% 
% % Find order of placefileds in selected category
% [~, sortvector] = max(cats{sortcategory}.spaceplot);
% [~, sortindices] = sort(sortvector);
% 
% for c = 1:length(data.metadata.categories)
%     if dFoT
%     FigHandle = figure('Position',[(c-1)*450,50,500,900]);
%     subplot(2,1,1);
%         imagesc(cats{c}.spaceplot(:,sortindices)');
%         title(['Category ' num2str(c)])
%         caxis([0.0 0.3]);
%     subplot(2,1,2);
%         imagesc(cats{c}.alltracks(:,sortindices)');
%         caxis([0.0 2.0]);
%     else
%     FigHandle = figure('Position',[(c-1)*450,50,500,900]);
%     imagesc(cats{c}.spaceplot(:,sortindices)');
%     title(['Category ' num2str(c)])
%     caxis([0.0 0.3]);
%     end
% end
% end
