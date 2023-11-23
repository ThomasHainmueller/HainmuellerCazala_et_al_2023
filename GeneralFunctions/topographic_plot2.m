function topographic_plot2(data,variable,varargin)
% Plot the extent of a selected variable in a given category color coded
% onto the outlines of the cell ROIs in a dataset. E.g. to visualize the
% magnitude of activity distribution in the local network.
% Arguments:
% variable: Any variable stored in the cell{} compartement
% plotcategory: set 1 if variable is independent of category.

args=struct('active',0,'SI',1.0,'PF',1.0,'region','',...
    'categories',1:length(data.metadata.categories),'combination','or',...
    'plotcategory',1,'range',[-1 1]); 
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

plotindices = findcells(data, args.plotcategory, args.active, args.SI,...
    args.PF, args.combination, args.region);

nplanes = data.cells{length(data.cells)}.plane; % Assume that last cell has highest plane number

% Plotting part
for p = 1:nplanes
    figure;
    hold on;
    for n = plotindices
        thisvar = getfield(data.cells{n},variable);
        thisvar = thisvar(args.plotcategory);
        thisroi=data.cells{n}.roi.mnCoordinates(:,:);
        
        % Create appropriate colour representation
        if thisvar <= args.range(1)
            thiscol = 0;
        elseif thisvar >= args.range(2)
            thiscol = 1;
        else
            thiscol = (abs(thisvar - args.range(1)))/(args.range(2)-args.range(1));
        end
        
        if isnan(thiscol) || data.cells{n}.plane~=p
            continue
        end
        
        thiscol = [1-thiscol,1-thiscol,1-thiscol];
        
        h=fill(thisroi(:,1),-thisroi(:,2),thiscol);
        set(h,'EdgeColor',[.75,.75,.75]);
    end
end   
end