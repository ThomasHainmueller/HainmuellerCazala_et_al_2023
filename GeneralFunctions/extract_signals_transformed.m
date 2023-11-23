function [data, mask] = extract_signals_transformed(video,rois,tform)
% Extract signals from rois in a video. Rois are given in a format returned
% by 'ReadImageJROI.m' from a FIJI.zip file and video should be a double or
% int16 3D array. Transform the ROIs by the tform given before signal
% extraction.

% Create binary roi masks
for r = length(rois):-1:1
    brois(:,:,r)=poly2mask(rois{r}.mnCoordinates(:,1),...
        rois{r}.mnCoordinates(:,2),size(video,1),size(video,2));
end

brois=imtransform(brois,tform,'XData',[1 size(video,2)],...
    'YData',[1 size(video,1)]);

% Loop over each frame of the video and extract all roi mean values
for f = size(video,3):-1:1
    plane=video(:,:,f);
    for r=size(brois,3):-1:1
        data(r,f)=mean(plane(brois(:,:,r)));
    end
end
mask=brois;
end