function [cleanROIs] = clearROIs(rois, minsize)
% Eliminate small rois from a logical array with multiple roi masks
if nargin < 2
    minsize = 10;
end

for n=1:size(rois,3)
    ROIsize = length(find(rois(:,:,n)));
    if ROIsize >= minsize
        roiselector(n)=true;
    end
end
cleanROIs = rois(:,:,roiselector);
end