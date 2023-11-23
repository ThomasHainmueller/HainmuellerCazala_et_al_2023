% Creates the average of the ROI in a g/r matrix (gorobj). Input arguments
% are the matrix, that is directly derived via "export matrix" from MES and
% the size of the ROI (an equaly sized region will be used for BG
% correction in the red channel). The function returns a one dimensional
% time series vector (=the G/R normalized average trace in the ROI).
% data = g/r matrix (export 2channel matrix from MES)
% ROIsize = size of ROI (in pixels) and background
% tRecording = total duration of the recording (in seconds)
% F0start, F0end = start and end of the normalization period (in seconds),
% important, it's a DELTAgreen over Red calculation
% fFilter = frequency (in Hz) to which results are downsampled (def 100Hz)
% to try this out, load green_red_matrix.mat and type t = roi_average(gr,4,0.25,1,4,100);

function trace = roi_average(data,ROIsize, F0start, F0end, tRecording, fFilter)

% Generate a one-dimensional vector as a mean of a time-series array data
[gr ROIind] = GoR(data,ROIsize,ROIsize,F0start,F0end,tRecording);
ROImatrix = gr(ROIind,:);
trace = mean(ROImatrix);

% Convert the raw trace vector to a timeseries, downsampled to fFilter
trace = resample(trace,fFilter*tRecording,length(trace));
trace = timeseries(trace);
trace.time = trace.time/fFilter;

end