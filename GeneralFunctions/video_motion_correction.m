function [T, refm] = video_motion_correction(video, nFramesRef, nIterations)
% Align a m x n x frames video using rigid alignment. First create a reference
% image from the initial f frames and then align all consecutive images by
% that. Return shifts in px for each frame.
if nargin < 3
    nIterations = 2; % number of frames used for 2nd alignment step (plus and minus)
end

if nargin < 2
    nFramesRef = 20;
end
[optimizer, metric] = imregconfig('monomodal');
optimizer.MaximumIterations = 1000;
optimizer.MaximumStepLength = .008;

% Create initial reference image for alignment
for n = 1:nFramesRef
    T{n} = imregtform(video(:,:,n),video(:,:,2),'translation',optimizer,metric);
    ref(:,:,n) = imwarp(video(:,:,n),T{n},'OutputView',imref2d(size(video(:,:,2))));
end

for it = 1:nIterations
    refm = mean(ref,3);
    
    parfor n = nFramesRef+1:size(video,3)
        T{n} = imregtform(video(:,:,n),refm,'translation',optimizer,metric);
        ref(:,:,n) = imwarp(video(:,:,n),T{n},'OutputView',imref2d(size(refm)));
    end
end

end