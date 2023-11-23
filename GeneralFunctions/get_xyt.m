function [green red metadata] = get_xyt(f,index)

% metadata:{frametimes xcords(framewise) ycords(framewise) xraw yraw}
%try
    %num_frames = get(f(index),1,'D3Size');
%catch
    num_frames = get(f(index),1,'FoldedFrameInfo');
    num_frames = num_frames.numFrames;
%end

% Get green data
for n = num_frames:-1:1
    green(:,:,n) = get(f(index),n,'IMAGE');
end

% Get red data
for n = num_frames:-1:1
    red(:,:,n) = get(f(index),num_frames+n,'IMAGE');
end
% blank emtpy frames
last_frame=find(all(green(:,end,:) == 0),1,'first')-1;

if isempty(last_frame)
    last_frame=num_frames;
end

green=green(:,:,1:last_frame);
red=red(:,:,1:last_frame);
% remove uneven background due to ambient light.
green = rm_ambientlight(green);
red = rm_ambientlight(red);

xcord_raw = gorconvert(get(f(index),1,'AUXin1Raw'));
ycord_raw = gorconvert(get(f(index),1,'AUXin2Raw'));

% Get the starttimes of the frames, and x- and y-coords for each frame.
start=get(f(index),1,'D3Origin');
interval=get(f(index),1,'D3Step');
frametimes=[start:interval:(last_frame-1)*interval+start];

for iter = length(frametimes):-1:1
    xstart = find(xcord_raw(:,1)>frametimes(iter),1,'first');
    xend = find(xcord_raw(:,1)>frametimes(iter)+interval,1,'first');
    xcord(iter)=mean(xcord_raw(xstart:xend,2));
    ycord(iter)=mean(ycord_raw(xstart:xend,2));
end

metadata{1}=frametimes;
metadata{2}=xcord;
metadata{3}=ycord;
        
end