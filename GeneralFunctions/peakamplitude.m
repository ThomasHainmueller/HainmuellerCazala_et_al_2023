% peakamplitude takes a g/r matrix and extracts the peak of amplitude after
% FOend in a ROI of specified width after normalizing to the red channel
% and filtering by fFilter. Variable description see roi_average.

function peak = peakamplitude(data,ROIsize, F0start, F0end, tRecording, fFilter)
tr = roi_average(data,ROIsize, F0start, F0end, tRecording, fFilter);
% cut of the start of the trace (usually large amplitude by shutter
tr = getsamples(tr,[F0end/tRecording*length(tr.time):length(tr.time)]);
peak = max(tr);
end