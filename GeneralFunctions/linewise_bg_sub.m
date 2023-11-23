function subtracted = linewise_bg_sub(array)
% Fraction of pixels per line assumed to contain background signal
bgfract = 0.9;
darkfract = 0.03;

% Find the darkest overall frames in the framescan.
[~,darkindices] = sort(mean(mean(array,1),2));
darkindices = darkindices(1:round(darkfract*length(darkindices)));
% Modified 160428; use only bgfract of darkest pixels to calculate line bg
% bg = mean(mean(array(:,:,darkindices),3),1);
bgimg = mean(array(:,:,darkindices),3); % background image (x,y)
% Get one global background value for each line reflecting the mean of the
% bgfract% darkest pixels
for h=1:size(array,2)
    bgline = sort(bgimg(:,h));
    bg(1,h)=mean(bgline(1:round(bgfract*length(bgline))));
end

for frame = 1:size(array,3)
    for h = 1:size(array,2)
        % get current line from frame
        line = array(:,h,frame);
        % get background value as [fraction] of the lowest intensity pixels
        line = sort(line);
        thisbaseline = mean(line(1:round(bgfract*length(line))));
        % subtract the basline change between the current and darkes frames
        % of the entire run.
        array(:,h,frame) = array(:,h,frame)-thisbaseline+bg(h);
    end
end
subtracted = array;
return
end

% Deprecated 160204
% offset = 1000;
% 
% for frame = 1:size(array,3)
%     for h = 1:size(array,2)
%         % get current line from frame
%         line = array(:,h,frame);
%         % get background value as [fraction] of the lowest intensity pixels
%         line = sort(line);
%         bg = mean(line(1:round(bgfract*length(line))));
%         % subtract background from line; offset is added as to correct for 
%         % zero setting of high noise values in the integer array.
%         array(:,h,frame) = array(:,h,frame)+offset-bg;
%     end
% end
% subtracted = array;
% return
