function signals = pullsignals(data, rois, tform)
% Pull signals for all specified rois for a given plane from an imaging
% dataset.
% data: m x n x frames 3D uint16
% rois: binary mask m x n x nrois 3D logical
% tform: Can be affine2d() or struct (maketform()) type 2D affine transform
% Loop over each frame of the video and extract all roi mean values

%-------------TODO: USE TINV FOR CALLING THIS FKT--------------

% Use gpuArray if enough memory is available on gpu

% g = gpuDevice();
% s = whos('data');
% 
% if 2*s.bytes <= g.AvailableMemory % Check enough memory is available on gpu
%     try
%         data = gpuArray(data);
%         gpuprocc = true;
%     catch
%         gpuprocc = false;
%     end
% else
%     gpuprocc = false;
% end
%
% MAJOR UPDATE (181026): Corrected bug in data extraction (use logical
% mask for indexing instead of m,n coordinates).

gpuprocc = false;

for r = length(rois):-1:1
    brois(:,:,r)=poly2mask(rois{r}.mnCoordinates(:,1),...
        rois{r}.mnCoordinates(:,2),size(data,1),size(data,2));
end

% Apply spatial transformation to rois
if isstruct(tform)
    tform = affine2d(tform.tdata.Tinv);
end
ref = imref2d(size(brois));
brois = imwarp(brois,tform,'OutputView',ref);

for r = size(brois,3):-1:1
    % Extract data frame-by-frame using logical indexing
    for f = 1:size(data,3)
        thisim = data(:,:,f);
        signals(f,r) = mean(thisim(brois(:,:,r)));
        clear thisim
    end 
    
%     thissignals = data(mask); %Alternate approach idea
%     thissignals = reshape(thissignals,[],size(data,3));
%     signals(:,r) = mean(thissignals,1);
%     clear thissignals mask

%     [m,n]=find(brois(:,:,r)); %Corrected 181026
%     signals(:,r) = mean(mean(data(m,n,:),2),1);
%    figure; C=imfuse(mean(data,3), brois(:,:,r)); imshow(C);
end

if gpuprocc
    signals=gather(signals);
end

end