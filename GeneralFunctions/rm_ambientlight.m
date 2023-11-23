function subtracted = rm_ambientlight(array)
% "Homogenizes" the background by replacing the mean of each line with the
% mean of this line in the %dark_periods.

% This is the fraction of time for which no backgroundlight is assumed
dark_fract = 0.1;

% Find darkest frames.
intensities=median(median(single(array),1),2);
intensities=reshape(intensities,size(intensities,3),1,1);
int_sort=sort(intensities);
cutoff=int_sort(round(0.1*length(intensities)));
darkest=intensities<=cutoff; %logical with the darkest frames
clear intensities int_sort cutoff;

test = array(:,:,darkest);
bg_lines = median(median(single(test),3),1);
%plot(reshape(bg_lines,length(bg_lines),1,1));

for frame = size(array,3):-1:1
    for h = size(array,2):-1:1
        % get current line from frame
        line = array(:,h,frame);
        % get median of current line
        bg = median(single(line),1);
        % subtract background from line; offset is added as to correct for 
        % zero setting of high noise values in the integer array.
        BGval = bg_lines(1,h,1)-bg;
        if BGval>0
            subtracted(:,h,frame) = array(:,h,frame)+uint16(BGval);
        else
            subtracted(:,h,frame) = array(:,h,frame)-uint16(-BGval);
        end
    end
end

end