function [stack]=imfilt_2D(stack)
% Use a 2D filtering methond, e.g. median or mean on an image stack.
tic
w=2; % window size

% Initialize and set up parameters
ymax=size(stack,1);
xmax=size(stack,2);
numFrames=size(stack,3);

% Apply filtering to each plane of the stack
for pl = 1:numFrames
    sprintf('Processing plane %d',pl);
    thispl = stack(:,:,pl);
    for y=1+w:ymax-w
        for x=1+w:xmax-w
            pix=reshape(thispl(y-w:y+w,x-w:x+w),1,(2*w+1)^2);
            stack(y,x,pl) = median(pix);
        end
    end
end
toc
end

function filt = imfilt_plane(img)
global w xmax ymax
% Apply filtering to individual plane
%filt=zeros(size(img)); % duplicate plane
1+w
ymax-w
xmax-w
for y=1+w:ymax-w
    for x=1+w:xmax-w
        pix=reshape(img(y-w:y+w,x-w:x+w),1,(2*w+1)^2);
        filt(y,x) = median(pix)
    end
end
end