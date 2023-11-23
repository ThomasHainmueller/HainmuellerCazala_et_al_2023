function topographic_plot(Cell,category,PVdata,PCdata,variable,range)
% Plot the Correlations between PVcells and Principal neurons of the same
% dataset on a topographic map. The function can plot any parameter on the
% topographic map, but is currently only implemented for this one purpose.
% variable: 'correlation'
if nargin < 5
    range=[-0.2,0.2];
end
if nargin < 4
    variable='correlation';
end

if strcmp(variable,'correlation')
    % Plot all PC rois in a color representing their correlation to the
    % selected PV cell.
    figure;
    hold on;
    for PCno=1:length(PCdata.cells)
        thiscorr=PVdata.cells{Cell}.categories{category}.PCcorr(PCno,1);
        thisroi=PCdata.cells{PCno}.roi.mnCoordinates(:,:);
        % Generate appropriate color for plotting this roi
        if thiscorr<=0
            col = thiscorr/range(1); % blue for negative values
            if col>1
                col=1; % Exclude invalid plotting values
            elseif isnan(col)
                col=0;
            end
            col = [1-col,1-col,1];
        else
            col = thiscorr/range(2); % red for positive values;
            if col>1
                col=1; % Exclude invalid plotting values
            elseif isnan(col)
                col=0;
            end
            col = [1,1-col,1-col];
        end
        h=fill(thisroi(:,1),thisroi(:,2),col);
        set(h,'EdgeColor',[.75,.75,.75]);
    end
    
    for PVno=1:length(PVdata.cells)
        thiscorr=PVdata.cells{Cell}.categories{category}.PVcorr(PVno,1);
        thisroi=PVdata.cells{PVno}.roi.mnCoordinates(:,:);
        % Generate appropriate color for plotting this roi
        if thiscorr<=0
            col = thiscorr/range(1); % blue for negative values
            if col>1
                col=1; % Exclude invalid plotting values
            elseif isnan(col)
                col=0;
            end
            col = [1-col,1-col,1];
        else
            col = thiscorr/range(2); % red for positive values;
            if col>1
                col=1; % Exclude invalid plotting values
            elseif isnan(col)
                col=0;
            end
            col = [1,1-col,1-col];
        end
        h=fill(thisroi(:,1),thisroi(:,2),col);
        set(h,'EdgeColor',[0,0,0]);
    end    
end