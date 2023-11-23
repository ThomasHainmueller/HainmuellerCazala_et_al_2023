function [green red metadata] = get_framescan(f,index,shift)
if nargin <3
    shift = 2;
end
% metadata:{frametimes xcords(framewise) ycords(framewise) xraw yraw}

% evalin('base',['hu =' num2str(23)]);
green_raw = get(f(index),1,'IMAGE');
red_raw = get(f(index),2,'IMAGE');
% gorconvert.m required for this.
try
    xcord_raw = gorconvert(get(f(index),1,'AUXin1Raw'));
    ycord_raw = gorconvert(get(f(index),1,'AUXin2Raw'));
catch
    xcord_raw = gorconvert(get(f(index),1,'AUXi1'));
    ycord_raw = gorconvert(get(f(index),1,'AUXi2'));
end

metadata{4}=xcord_raw;
metadata{5}=ycord_raw;

params = get(f(index),1,'FoldedFrameInfo');
nframes = params.numFrames - 1;
% last frame is usually not complete - omit!
begin = params.firstFramePos;
heigth = size(green_raw,1);
width = params.numFrameLines;

green=[];
red=[];

frametimes(1) = params.firstFrameStartTime;
for iter = 1:nframes-1
    gframe = green_raw(:,begin+width*(iter-1):begin+width*iter-1);
    rframe = red_raw(:,begin+width*(iter-1):begin+width*iter-1);
    
    % evaluate whether the last scanline in this frame was still acquired
    if gframe(end,:)~=0
        green(:,:,iter)=gframe;
        red(:,:,iter)=rframe;
        frametimes(iter+1)=frametimes(1)+iter*params.frameTimeLength;
        % Logical mask to get the x-y coords during the current frame
        x_indices =(xcord_raw(:,1)>frametimes(iter))&(xcord_raw(:,1)<frametimes(iter+1));
        y_indices =(ycord_raw(:,1)>frametimes(iter))&(ycord_raw(:,1)<frametimes(iter+1));
        % Get mean x-y coordinate value during current frame
        xcord(iter)=mean(xcord_raw(x_indices,2));
        ycord(iter)=mean(ycord_raw(y_indices,2));
    else
        break
    end  
end
% corrects for bidirectional scanning error, requires shiftlines.m
green = shiftlines(green,shift);
red = shiftlines(red,shift);
metadata{1}=frametimes;
metadata{2}=xcord;
metadata{3}=ycord;
        
xrawrange = (xcord_raw(:,1)<frametimes(iter));
metadata{4}= xcord_raw(xrawrange,:);
        
yrawrange = (ycord_raw(:,1)<frametimes(iter));
metadata{5}= ycord_raw(yrawrange,:);

%green = green(:,begin:begin+nframes*width-1);
%red = red(:,begin:begin+nframes*width-1);

% corrects for bidirectional scanning error, requires shiftlines.m
%green = shiftlines(green,shift);
%red = shiftlines(red,shift);

%green = permute(green,[2 1 3]);
%red = permute(red, [2 1 3]);
%return
end