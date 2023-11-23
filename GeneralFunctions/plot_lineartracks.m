function plot_lineartracks(dataset,cellno,normalized)
% Plots lineartrack acquired activity grouped to category. Note that the
% number of cells is 3 smaller than the size of the dataset as 1-3 hold
% frametimes, x- and y coordinates.

if nargin < 3
    normalized = true;
end

res = {};
ft = {};

% Append all relevant Data to cell arrays
for n = length(dataset):-1:1
    thiscat = dataset{n}.category;
 
    %thisy = thisy/max(max(thisy));
    
    res{thiscat,end+1}=dataset{n}.data(3,:);
    %res{thiscat,end+1}=dataset{n}.data(cellno,:)-1.0;
    if normalized
        res{thiscat,end+1} = normtobaseline(dataset{n}.data(cellno,:),0.1);
    else
        res{thiscat,end+1}=dataset{n}.data(cellno,:);
    end
    ft{thiscat,end+1}=dataset{n}.data(1,:);
    ft{thiscat,end+1}=dataset{n}.data(1,:);
end

% Convert cellarrays to plottable matlab arrays
for m = size(res,1):-1:1
     % Cellarray from res is converted to a double variable which has the
     % legth of the longest trace in this category. Shorter traces are
     % padded with NaN.
     thisres=catuneven(res(m,:),NaN);
     thisft = catuneven(ft(m,:),NaN);
        
     thisres(~any(~isnan(thisres),2),:)=[];
     thisft(~any(~isnan(thisft),2),:)=[];
        
     pres{m} = thisres;
     pft{m} = thisft;
     clear thisres thisft
end

dFoY = {};
% Make dF/ycoordinate plots for each trace ordered by category
for m = size(res,1):-1:1
    % bins describes the discrete bins, by which the y coordinates are
    % grouped. Cave, y values in res are normalized to 1.0!
    bins = [0.05:0.05:2.1];
    thisdFoY = double.empty(size(pres{m},1)/2,size(bins,2)-1,0);
    for o = size(pres{m},1)/2:-1:1
        thisdFoY(o,:,1) = SBdiscretize(pres{m}(2*o,:),pres{m}(2*o-1,:),bins);
    end  
    dFoY{m}=thisdFoY;
    clear thisdFoY
end

% Create the final plot, defaults to 6x4 subplots (12 categories, traces
% and image each.

for subind = size(pres,2):-1:1
    if size(pres{subind},1) > 1
        figure
        subplot(2,2,[1 3])
        % All traces on top of each other, each as y coord / trace pair
            for iter = size(pres{subind},1):-1:1
                hold on;
                plot(pft{subind}(iter,:),(pres{subind}(iter,:)-iter));
            end
            %plot3(pft{subind}',pylines{subind}',pres{subind}')
            title(['Category ' num2str(subind)])
        subplot(2,2,2)
        % dF over y coord for each run in this category
            imagesc(dFoY{subind})
            % TODO: reset to -0.5 1.0
            caxis([-0.1 1.0])
        subplot(2,2,4)
        % Mean +/- std for dF over y coord for all runs together
            plot(nanmean(dFoY{subind}),'k')
            hold on
            plot(nanmean(dFoY{subind})-nanstd(dFoY{subind}),'r')
            hold on
            plot(nanmean(dFoY{subind})+nanstd(dFoY{subind}),'r')

    else
        warning(['category ' num2str(subind) ' empty'])
    end
end

