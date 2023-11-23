function [map] = map_data(trace, xcor, ycor)
% Create an array of data over space. Like 'spatialplot' but returns.
ycor = 10*ycor;
xcor = 10*xcor;
ar = {};

% Cellarray holding all values that were acquired in one xy-position.
for iter=[1:length(xcor)]
    try
        ar{round(xcor(iter)),round(ycor(iter))} = cat(1,thispx,trace(iter));
    catch
        ar{round(xcor(iter)),round(ycor(iter))} = trace(iter);
    end
end

%Create a heat mapped image of the mean Calcium signal in each xy-pos.
map = zeros(round(max(xcor)),round(max(ycor)));
for h = [1:size(map,1)]
    for w = [1:size(map,2)]
        map(h,w) = mean(ar{h,w});
    end
end
end