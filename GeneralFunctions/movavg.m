function smooth = movavg(trace, windowsize)
% Create a smoothed trace using a moving average procedure. Ends are
% averaged only with half the window size.

for n = length(trace):-1:1
    if n-windowsize/2 < 1
        window = 1:round(n+windowsize/2);
    elseif n+windowsize/2 > length(trace)
        window = round(n-windowsize/2):length(trace);
    else
        window = round(n-windowsize/2):round(n+windowsize/2);
    end
    smooth(n) = mean(trace(window));
    clear window;
end

end