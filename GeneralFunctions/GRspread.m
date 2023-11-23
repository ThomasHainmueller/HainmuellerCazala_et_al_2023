% GRspread takes a GoR normalized matrix and extracts the peak amplitude
% from each pixel along the scan line allowing the user to construct a
% histogram of the signal spread. It filters the data by f Filter and
% performs a spatial average (boxcar) with a windowsize defined by +/- SpAvg.
% This function is usefull only for LONGITUDINAL LINESCAN measurements.
% Cave: ends of the line are clipped by SpAvg

function peakamplitudes = GRspread(data, F0start, tRecording, fFilter, SpAvg)
peakamplitudes = [];
for i = 1+SpAvg:size(data,1)-SpAvg
    trs = data([i-SpAvg:i+SpAvg],:);
    tr = mean(trs);
    % "filter function" from roi average
    tr = resample(tr,fFilter*tRecording,length(tr));
    % cut of the shutter movement, extract peak value for each pixel (tr)
    tr = tr(F0start*fFilter:tRecording*fFilter);  
    %peak = max(tr);
    peakamplitudes = [peakamplitudes max(tr)];
    % debugging
    plot(tr);
    hold on;
end
end