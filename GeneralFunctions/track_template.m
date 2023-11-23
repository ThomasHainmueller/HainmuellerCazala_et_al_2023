function [data, outvid] = track_template(video, template, varargin)
% An object on a square image is tracked in each frame of a video file and
% the x/y translation plus rotation are returned in a 3 x nframes matrix.
% Additionally, the result can be displayed as an additive video, where the
% matched template is superimposed onto each frame of the output video. The
% template should contain a *square* top view of the circular arena with the 
% centre aligned with the image. For adequate interpretation of the mouse
% loaction in the MHC, the video should be aligned so that the mouse is in
% it's centre
% 
% Designed to identify the position of the mobile homecage in a birdseye
% view video. Output: 7 x nframes matrix, representing x,y and rotation of
% the MHC, x,y position of the mouse relative to the MHC centre plus a
% vectorial representation (angle, distance) of the same.
%
% Recommended preprocessing: Make image of MHC in the middle of the FOV,
% crop MHC out, Crop video to the relevant area, but make sure that the
% entire MHC is visible at all times!

args=struct('folder',[]);

% Overwrite default parameters if required
for pair = reshape(varargin,2,[])
    if isfield(args,pair{1})
        args.(pair{1})=pair{2};
    else
        error('Input argument does not exist\n')
    end
end

% Preprocessing steps: Inversion and bg subtraction
template = iminvert(template);
%bgim = uint8(mean(video,3));
for f = size(video,3):-1:1
%    video(:,:,f) = iminvert(video(:,:,f))-bgim;
    video(:,:,f) = iminvert(video(:,:,f));
end

% Create circular mask to avoid corner problem when rotating
[xx, yy] = meshgrid(1:size(template,2), 1:size(template,1));
mask = (yy - round(size(template,1)/2)).^2 ...
    +(xx - round(size(template,2)/2)).^2 <= round(size(template,1)/2).^2;

% Initialize output video (RGB)
outvid(:,:,1,:) = reshape(video,[size(video,1),size(video,2),1,size(video,3)]);
outvid(:,:,2:3,:) = 0;

for f = size(video,3):-1:1
    % Identify the most likely x/y position of the template
    if exist('thistemp','var')
        % Use previous (rotated) template if available
        c = normxcorr2(thistemp, video(:,:,f));
    else
        c = normxcorr2(template, video(:,:,f));
    end

    % Implicitly assume that the disk is completely in the FOV at all times
    c([1:size(template,2)+1, end-size(template,2):end],:) = 0;
    c(:,[1:size(template,1)+1, end-size(template,1):end]) = 0;
    
    [~, imax]= max(c(:));  
    [ypeak, xpeak] = ind2sub(size(c), imax(1));
    offset = [(xpeak-size(template,2)), (ypeak-size(template,1))];
    
    % Get area containing the template from this frame
    match = video(round(offset(2)+1):round(offset(2)+size(template,2)),...
        round(offset(1)+1):round(offset(1)+size(template,1)),f);
    match = match.*uint8(mask);
    
    % Find best rotation 1° increments
    for angle = 1:360
        thistemp = center_rotate(template, angle); 
        thistemp = thistemp.*uint8(mask);
        rotcorr(angle) = corr2(thistemp, match);
    end
   
    % Create optimally rotated template
    [~, bestrot] = max(rotcorr);
    thistemp = center_rotate(template, bestrot);

    % Identify optimal x/y shift using rotated template
    c = normxcorr2(thistemp, video(:,:,f));
    
    c([1:size(template,2)+1, end-size(template,2):end],:) = 0;
    c(:,[1:size(template,1)+1, end-size(template,1):end]) = 0;
    
    [~, imax]= max(c(:));  
    [ypeak, xpeak] = ind2sub(size(c), imax(1));
    offset = [(xpeak-size(template,2)), (ypeak-size(template,1))];   
    
    % Create allocentric coordinates of the mouse (assumed as video centre)
    % relative to the MHC
    
    % Distances between video centre and MHC centre
    xdist = offset(1) + size(template,2)/2 - size(video,2)/2;
    ydist = offset(2) + size(template,1)/2 - size(video,1)/2;
    absdist = sqrt(xdist^2 + ydist^2);
    
    % Angle between template 'north' and video centre (RB)
    RB = mod(atan(ydist/xdist)/(2*pi) *360 +180 -bestrot, 360);
    
    % Mouse position relative to MHC centre
    mousex = cosd(RB)*absdist;
    mousey = sind(RB)*absdist;
    
    % Store values for output
    data(1:7,f) = [xdist, ydist, bestrot, mousex, mousey, RB, absdist];   
    %data(1:7,f) = [offset(1), offset(2), bestrot, mousex, mousey, RB, absdist];   

    % Add matched template to the 2nd channel of output video
    outvid(offset(2)+1:offset(2)+size(template,2),...
        offset(1)+1:offset(1)+size(template,1),2,f) = thistemp;

end
end

function rotated = center_rotate(image, angle)
% Rotate an image around its center and crop to the original image size.
rotated = imrotate(image, angle);
% Crop and mask
xoffset = floor((size(rotated,2)-size(image,2))/2);
yoffset = floor((size(rotated,1)-size(image,1))/2);
rotated = rotated(yoffset+1:yoffset+size(image,1),...
    xoffset+1:xoffset+size(image,2));
end

function inv = iminvert(im)
% Invert a uint8 image
inv = uint8(ones(size(im))*255)-im;
end