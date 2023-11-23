function [registered, tform] = xcorr2_reg_median(fixedstack, moving, maxshift)
% Align a moving image to a stack of reference images. Determine the best
% fit for each individual image and return the median transformation matrix
% and aligned image as the best fit to all images in the stack.

for n = size(fixedstack,3):-1:1
    % Register moving with each plane of refstack, infer x and y shifts
    [~,tform]=xcorr2_reg(fixedstack(:,:,n),moving,maxshift);
    xshift(n)=tform.tdata.T(3,1);
    yshift(n)=tform.tdata.T(3,2);
end

xshift=median(xshift);
yshift=median(yshift);

% Make new tform with the median shift values
Atform=[1 0 0; 0 1 0; xshift yshift 1];
tform = maketform('affine',Atform);

% Apply image transformation
registered = imtransform(moving,tform,'XData',[1 size(fixedstack,2)],...
    'YData',[1 size(fixedstack,1)]);
end