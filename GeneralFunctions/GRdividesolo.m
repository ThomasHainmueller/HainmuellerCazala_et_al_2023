% Splits a long g/r array to single sweeps of defined length. Returns a
% cell with the individal g/r arrays for each sweep.
% data = g/r array (derived from MES->export matrix
% tRecording = total recording duration of data (s)
% interval = length of one sweep (s) => try 3.256
% delay = delay from start of recording to first sweep (s)
% Soloversion is takes one solitary matrix as input. For cells use GRdivide


function GRcell = GRdividesolo(data, tRecording, interval, delay)

if nargin < 4
    delay = 0.75;
end 

if nargin < 3
    interval = 3.23182;
end

if nargin < 2
    tRecording = 40;
end

% calculate interval and delay as indices (adjust to sampling rate); note
% that it is inherently assumed that the time domain is the longest axis of
% the data array!
delay = round(delay/tRecording*length(data));
interval = round(interval/tRecording*length(data));

% n is the number of sweeps into which data is divided
n = fix((length(data)-delay)/interval);

for c = 1:n
    first = delay+(c-1)*interval;
    last = delay+c*interval;
    GRcell{c} = data(:,first:last,:);
end
end
