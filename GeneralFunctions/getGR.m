% Extract the G/R matrices from subfiles in a MESfile (mestaghandle). Get
% the file via f=mestaghandle('isf'); Either create one file with all
% traces that you want to analyze or parse a vector containing the
% subfileindices to getGR. Returns a cell containing all the G/R matrices.

% Tipp: To merge cells from different files, use data = {cell1{:} cell2{:}}

function GRcell = getGR(MESfile, traces)

if nargin < 2
    % get all G/R matrices in this file
    traces = [1:length(MESfile)];
end 

%GRcell = cell(1,length(traces));

for c = [1:length(traces)]
    % extract Matrices from file, c is the subfile index
    index = traces(c);
    GR = double.empty();
    GR(:,:,1) = get(MESfile(index), 1, 'IMAGE');
    GR(:,:,2) = get(MESfile(index), 2, 'IMAGE');
    GRcell{c} = GR;
end
end