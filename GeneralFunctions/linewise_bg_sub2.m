function array = linewise_bg_sub2(array, fbg)
% Conceptual idea: Use the 'fbg' darkest pixels of each line's timeaverage
% and determine their mean timecourse as a dynamic background signal.
% subtract this dynamic background for each pixel in the corresponding line
% and set the baseline value to the mean of this dynamic signal (this to
% avoid distorsions by <0 values in noisy signals).

if nargin < 2
    fbg = .1; %Fraction of pixels per line which contain just background
end

bgimg = uint16(mean(array,3));

for l = 1:size(array,1)
    % Iterate over lines, find bg signal for each line
    [~,bgind] = sort(bgimg(l,:));
    bgind = bgind(1:round(fbg*length(bgind)));
    bgsig = uint16(mean(array(l,bgind,:),2));
    
    % Iterate over pixels in each line, subtract dynamic bg and replace w.
    % offset value (= mean of dynamic bg).
    for p = 1:size(array,2)
        array(l,p,:) = array(l,p,:) + nanmean(bgsig,3) - bgsig;
    end
end  
end

