function Ispatial = spatial_info(trace, y, bins, framerate)
% Calculates the amount of spatial information contained in a trace (e.g. a
% Calcium signal) about the variable y (usually linear track distance).
% Feed this with the concatenated 'motion only' periods of the y and
% calcium transients. Also, use the 'significant-extracted' calcium
% transients!
%
% Ref. Skaggs et al., 1993. 
if nargin<4
    framerate=1; % Return bits/frame
end

if nargin<3
    bins = 0.1:0.05:2.1; % For Lineartrack4A and descendants.
end

% Get mean dF/F (approximates overall mean firing rate).
% Set values smaller 0 to 0(in this case this are very rare outliers in
% identified (positive) calcium transients.
trace(trace<0)=0;

% Introduced 3/25/2023: Remove values outside of bin range prior to
% calculation. These are differently handled by 'hist' and 'discretize' and
% will lead to erroneous results otherwise!
% 6/10/2023: Fixed bug here that shifts y and trace!!!
goodidcs = y>bins(1) & y<= bins(end);
y = y(goodidcs);
trace = trace(goodidcs);

meanrate = mean(trace);
step=(bins(end)-bins(1))/(length(bins)-1);
yhist = hist(y,bins(1)+step/2:step:bins(end)-step/2);
trhist = SBdiscretize(trace,y,bins);

for x = length(trhist):-1:1
%     thisp = yhist(x); 3/25/2023
    thisp = yhist(x)/length(trace);
    thisrate = trhist(x);
    
    % Spatial information, Skaggs et al., 1993
    Ispatial(x) = thisrate*log2(thisrate/meanrate)*thisp;
end
Ispatial = nansum(Ispatial)*framerate;
%Ispatial = nansum(Ispatial)/length(trace)*framerate; 3/25/2023
end