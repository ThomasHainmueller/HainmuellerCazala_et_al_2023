function dataset = add_planes(dataset,varargin)
% Take existing dataset and add the plane in which each cell was recorded.
% Comes automatically with datasets created after 170913.

args=struct('folder','','roiname','manual'); 
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

% Load the ROIs; Assume dataset has as many planes as roifiles exist
planeindex = [];
ROIfiles = ls(['*',args.roiname,'*']);
for f = 1:size(ROIfiles,1)
    rois{f} = ReadImageJROI(sprintf('%s\\%s_plane%d.zip',cd(),args.roiname,f));
    %Nrois(f) = length(rois{f});
    planeindex(end+1:end+length(rois{f})) = f;
end

if ~length(planeindex) == length(dataset.cells)
    warning('Provided ROI set does not match with dataset');
    return;
else
    for n=1:length(planeindex)
        dataset.cells{n}.plane = planeindex(n);
    end
end

end