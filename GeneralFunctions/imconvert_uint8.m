 function [uint8stack, varargout] = imconvert_uint8(stack,min, max)
% Convert an image stack to uint8. Min and max give the lower and upper
% boundary (in % of the total intensity distribution) of the range that is
% converted. Adjust this in order to exclude background and saturated
% pixels.

if nargin < 3
    max = 0.999;
end
if nargin < 2
    min = 0.001;
end

hist = reshape(stack,(size(stack,1)*size(stack,2)*size(stack,3)),1,1);
hist = sort(hist);
minval=hist(round(min*length(hist)));
maxval=hist(round(max*length(hist)));
clear hist;

stack(stack<=minval)=minval;
stack(stack>=maxval)=maxval;

uint8stack=uint8(stack/(maxval/255));
varargout{1}=minval;
varargout{2}=maxval;

end