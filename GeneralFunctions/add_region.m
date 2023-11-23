function dataset = add_region(dataset,varargin)
% Needs an exisiting dataset plus a ROI for every plane that defines the
% region where cells should be located, e.g. 'SupCA1_plane1',... This ROI
% should exist for every plane of the dataset but must not contain cells. 

args=struct('folder','','roiname',''); 
% Make sure to mention the ROIs according to which the dataset has been
% created!

% Overwrite default parameters if required
for pair = reshape(varargin,2,[])
    %args.(pair{1}) = pair{2};
    if isfield(args,pair{1})
        args.(pair{1})=pair{2};
    else
        error('Input argument does not exist\n')
    end
end

if args.folder
    cd(folder);
end

ROIfiles = ls(['*',args.roiname,'*']);
for f = 1:size(ROIfiles,1)
    % Get ROI for the region selected
    thisroi = ReadImageJROI(sprintf('%s\\%s_plane%d.roi',cd(),args.roiname,f));
    thisroi = thisroi.mnCoordinates;
    
    for n = 1:length(dataset.cells)
        if dataset.cells{n}.image_plane == f
            center = mean(dataset.cells{n}.roi.mnCoordinates,1); %Find the cells central point
            
            % If cell is in region, add label 'region'
            if inpolygon(center(1),center(2),thisroi(:,1),thisroi(:,2))
                dataset.cells{n}.region = args.roiname; %Cave, will overwrite
            end
        end
    end
end

end