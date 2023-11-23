function nonoverlap = RemoveOverlapping(rois, maxoverlap)
if nargin < 2
    maxoverlap = 0.8;
end

n=1;
while n <= size(rois, 3)
    % size of the current ROI
    RoiSize = length(find(rois(:,:,n)));
    % Loop over all remaining ROIs to check overlap
    for m = 1:size(rois,3)
        if m ~= n
            % Number of pixels shared between the two ROIs
            %shared = rois(:,:,n).*rois(:,:,m);
            %shared = shared > 0;
            %overlap = length(find(shared))/RoiSize;
            overlap = length(find(rois(:,:,n) & rois(:,:,m)))/RoiSize;
            if overlap > maxoverlap
                % Remove the currently tested ROI from set
                rois=rois(:,:,[1:n-1 n+1:size(rois,3)]);
                n = n-1; % set back counter, otherwise a roi is skipped
                break
            end
        end
    end
    n=n+1;
end
nonoverlap = rois;  
end