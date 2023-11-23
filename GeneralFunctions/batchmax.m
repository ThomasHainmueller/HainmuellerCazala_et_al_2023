% batchmax takes a cell of g/r matrices (directly from MES) and returns the
% maximum values of the dG/R normalized and filtered ROIs as a vector. For
% argument descriptions see roi_average/GoR respectively. Use batchmax as
% follows: extract all g/r matrices of one experiment from MES using the
% 'export matrix' command, create 'datacell = {gr1 gr2 gr3 ... grn}' and
% enter datacell to batchmax.

% CAVE: The automatic ROI selection works only if the brightest pixels in
% the red channel are the ROI (typically the case only in TRANSVERSE
% LINESCANS.

% For GRdivide splitted traces the following works:
% res = batchmax(data,4,0.001,0.4,3.2,100);

function [bmax traces] = batchmax(datacell,ROIsize, F0start, F0end, tRecording, fFilter)
bmax = [];
traces = [];
for i = 1:length(datacell)
    % appending maxvalues for each g/r matrix ROI to bmax.
    p = peakamplitude(datacell{i},ROIsize, F0start, F0end, tRecording, fFilter);
    bmax = [bmax p];
    % optional plotting of the ROI traces for control.
    thistrace = roi_average(datacell{i},ROIsize, F0start, F0end, tRecording, fFilter);
    traces = [traces thistrace];
    plot(thistrace);
    hold on;
end
bmax = transpose(bmax);
end