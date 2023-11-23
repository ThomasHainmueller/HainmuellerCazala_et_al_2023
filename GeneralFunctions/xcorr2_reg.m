function [registered, tform] = xcorr2_reg(fixed, moving, maxshift, rotation)
% Uses xcorr2 (2D cross correlation) to infer the translation of two
% images (fixed and moving) necessary to achieve the highest corralation
% between the two images (registration). Returns the registered image as
% well as the transformation Matrix (tform) that can be applied hereafter
% to other images (e.g. second channel or image stacks). Useful to register
% a signal channel (e.g. GCaMP) depending on the translation of a
% stationary (e.g. TdT) background image. In addition, a range of values
% for rotation angles can be specified to test for the angle that yields
% the maximum overlap. Angles ar give as degrees°.

% CAVE: This may not work with different sized images!
% Subtract mean
registered=moving;
moving=moving-mean(moving(:));
fixed=fixed-mean(fixed(:));

if nargin < 4
    rotation = 0;
end

% 1. Find the values for translation and rotation that yield the largest
% cross correlation between the aligned images.

% Sequential calculation of xcorr for all given roation angles
for r=length(rotation):-1:1
    % Rotate image
    thisR=rotation(r)/180*pi;
    ROTtform = [cos(thisR) sin(thisR) 0; -sin(thisR) cos(thisR) 0; 0 0 1];
    ROTtform = maketform('affine',ROTtform);
    thismov = imtransform(moving,ROTtform);
    
    % Calculate cross correlation
    xcorim = xcorr2(fixed, thismov);
    centerx = ceil(size(xcorim,1)/2);
    centery = ceil(size(xcorim,2)/2);
    
    % Store value of the maximum achieved cross correlation
    maxcorr(r)=max(max(xcorim(centerx-maxshift(1):centerx+maxshift(1),...
        centery-maxshift(2):centery+maxshift(2))));
end

% Find the rotation value that yielded the largest xcorr value
bestrot=rotation(find(maxcorr==max(maxcorr(:))));

% Rotate moving with the best angle
maxR=bestrot/180*pi;
ROTtform = [cos(maxR) sin(maxR) 0; -sin(maxR) cos(maxR) 0; 0 0 1];
ROTtform = maketform('affine',ROTtform);
thismov = imtransform(moving,ROTtform,'XData',[1 size(fixed,2)],'YData',[1 size(fixed,1)]);

% find coordinates of the maximum spatial cross correlation, within the
% limits of the maximal allowed shift.

xcorim=xcorr2(fixed,thismov);
centerx = ceil(size(xcorim,1)/2);
centery = ceil(size(xcorim,2)/2);
mask = zeros(size(xcorim,1),size(xcorim,2));
mask(centerx-maxshift:centerx+maxshift,centery-maxshift:centery+maxshift)=1;
xcorim = xcorim.*mask;
[maxX,maxY]=find(xcorim==max(xcorim(:)));


% Calculate X and Y shift. 
xshift = maxX-size(thismov,1);
yshift = maxY-size(thismov,2);

% Create translation tform based on x and y shift and rotation
Atform = [cos(maxR) sin(maxR) 0; -sin(maxR) cos(maxR) 0; yshift xshift 1];
tform = maketform('affine',Atform);

% Apply tform to moving. Result has the same dimensions as fixed.
registered = imtransform(registered,tform,'XData',[1 size(fixed,2)],...
    'YData',[1 size(fixed,1)]);
end