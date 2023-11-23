function ROIs = FindRois(green)
% Find ROIs of highly correlated adjacent pixels of a given size which
% likely correspond to individual cells in a calcium-imaging raster scan.

% Minimum correlation to neighbours for beein considered at all
threshold = 0.25;
% Minimum size (pixels) for a ROI
min_size = 15;
% Maximum overlap between individual ROIs
max_overlap = 0.5;

% Get the local cross-correlation map
green = single(green);
ccim = CrossCorrImage(single(green));
green=uint8(green);
myfilter = fspecial('gaussian',[5 5], 0.5);
ccimfilt = imfilter(ccim, myfilter, 'replicate');
thresholdmask = ccimfilt > threshold;
ccimfilt = ccimfilt.*thresholdmask;

% Find regional maxima in the cross-correlation map.
centers = imregionalmax(ccimfilt,8);

% ROIs of neighbouring pixels which are highly correlated to their center
ROIs = CCRois(centers,green);
clear green;

% Remove ROIs which dont match min_size and max_overlap criteria
ROIs = clearROIs(ROIs,min_size);
ROIs = RemoveOverlapping(ROIs,max_overlap);
figure, imshowpair(max(ROIs,[],3),ccimfilt);

end