function data=append_ROIs(data,roipath)
% Load FIJI ROIs from a .zip file and append them to a dataset.

rois = ReadImageJROI(roipath);

for n = 1:length(rois)
    % Add ROI centers
    rois{n}.center=[round(mean(rois{n}.mnCoordinates(:,1))),...
        round(mean(rois{n}.mnCoordinates(:,2)))];
    % Append ROIs to dataset
    data.cells{n}.roi = rois{n};
end

end