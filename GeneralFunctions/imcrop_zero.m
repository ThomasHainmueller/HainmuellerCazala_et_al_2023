function [xcords,ycords] = imcrop_zero(stack)
% Crops rows and columns that are full with zeros in one image of a stack.
% Good i.e. for finding the common part of a series of aligned images. This
% function returns the coordinates for cropping instead the cropped image.

% Make logical mask with the rows and columns containing nonzero values in
% all images.
logicx = any(all(stack==0,2),3);
logicy = any(all(stack==0,1),3);

% Find start and stop x,y coordinates
[x,trash]=find(~logicx);
[trash,y]=find(~logicy);

% crop the original stack, first (but not last) x and y coordinate must be
% incremented by one.
cropped=stack(x(2:end),y(2:end),:);

end