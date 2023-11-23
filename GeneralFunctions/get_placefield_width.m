function data = get_placefield_width(data,plotting)
% Append width of a cells placefield to dataset. If multiple placefields
% exist, width will append the cumulative sum of the field widths.
% Field width is given as number of spatial bins.

if nargin<2
    plotting = false;
end

if plotting
    figure; hold on;
end

for n=1:length(data.cells)
    for c=1:length(data.cells{n}.categories)
        data.cells{n}.placefield_width(c)=...
            length(find(data.cells{n}.categories{c}.placefield));
        
        if plotting
            cval = c/(2*length(data.cells{n}.categories));
            plot(n,data.cells{n}.placefield_width(c),'o','color',[cval cval cval]);
        end
    end
end

end