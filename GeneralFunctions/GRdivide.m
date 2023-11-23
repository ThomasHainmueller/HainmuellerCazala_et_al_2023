% Takes a cell containing G/R arrays and splits all arrays within by the
% parameters specified. It returns a cell containing all the splitted
% G/R arrays in a row.
% data = g/r array (derived from MES->export matrix
% tRecording = total recording duration of data (s)
% interval = length of one sweep (s) => try 3.256
% delay = delay from start of recording to first sweep (s)

function GRcell = GRdivide(datacell, tRecording, interval, delay)

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

for a = [1:length(datacell)]
    data = datacell{a};
    thisdelay = round(delay/tRecording*length(data));
    thisinterval = round(interval/tRecording*length(data));

    % n is the number of sweeps into which data is divided
    n = fix((length(data)-thisdelay)/thisinterval);
    
    for c = 1:n
        first = thisdelay+(c-1)*thisinterval;
        last = thisdelay+c*thisinterval;
        
        % the GRcell will contain a total of length(datacell)*n entries.
        % This only works if each entry in datacell is divided in an equal
        % number of segments (as tRecording and interval are the same for
        % all segments, this should usually be the case).
        cellindex = (a-1)*n+c;
        GRcell{cellindex} = data(:,first:last,:);
    end
end
end
