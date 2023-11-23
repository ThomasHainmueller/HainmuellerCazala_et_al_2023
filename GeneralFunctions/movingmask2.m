function mask = movingmask2(trace, framerate, minspeed)
% Creates a mask for a given trace that is 1 if the value of the smoothed
% trace exceeds minspeed. Minspeed in m/s assuming that 1.0 Units of
% distance in 'trace' equals 2 meters of real-world distance (as in the
% ball-VR system).
% DEPENDENCIES: running_mean()

thresh = minspeed/(2*framerate); % 1.0 Units = 2 m real dist.
meanwindow = 5;

dtrace = sqrt(diff(trace).^2);
%dtrace = running_mean(dtrace, meanwindow);
dtrace = conv(dtrace, ones(meanwindow,1)/meanwindow, 'same'); %Running mean filter
mask = [0, dtrace'>thresh];
end