function [green red] = get_zstack(f,index)

% metadata:{frametimes xcords(framewise) ycords(framewise) xraw yraw}
num_frames = get(f(index),1,'D3Size');

% Get green data
for n = num_frames:-1:1
    green(:,:,n) = get(f(index),2*n-1,'IMAGE');
end

% Get red data
for n = num_frames:-1:1
    red(:,:,n) = get(f(index),2*n,'IMAGE');
end
% blank emtpy frames
last_frame=find(all(green(:,end,:) == 0),1,'first')-1;

if isempty(last_frame)
    last_frame=num_frames;
end

green=green(:,:,1:last_frame);
red=red(:,:,1:last_frame);

        
end