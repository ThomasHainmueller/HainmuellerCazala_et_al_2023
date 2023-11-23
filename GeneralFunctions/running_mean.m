% Running mean of a matrix m with windowsize w, rowwise, ends are averaged
% to w first/last entries.
function smooth = running_mean(m, w)
smooth = [];
last = size(m,2);

% loop through lines
for l = [1:size(m,1)]
    % treating ends
    smooth(l,[1:w+1]) = mean(m(l,[1:w]));
    smooth(l,[last-w:last]) = mean(m(l,[last-w:last]));
    % loop along line l
    for c = [w+1:last-w]
        smooth(l,c) = mean(m(l,[c-w:c+w]));
    end
end
end