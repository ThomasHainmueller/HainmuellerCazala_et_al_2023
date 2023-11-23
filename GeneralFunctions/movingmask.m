function mask = movingmask(trace, minspeed)
% Creates a mask for a given trace that is 1 if the value of the smoothed
% trace exceeds minspeed.

meanwindow = 2;

dtrace = sqrt(diff(trace).^2);
dtrace = running_mean(dtrace,meanwindow);
mask = [0, dtrace>minspeed];
end