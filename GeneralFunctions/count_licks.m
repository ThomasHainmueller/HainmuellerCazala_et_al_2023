function lickcount = count_licks(rawtrace,fs,frames)
% Use a raw analog recording trace from arduino capacitive sensing, extract
% the lick-events and return the rate for each frame. 'frames' is a
% steadily growing vector where x is the datapoint and y is the frame in
% which it was recorded.

rawtrace = remove_50Hz(rawtrace,fs); % Filter 50 HZ noise
[b,a] = butter(4,.05,'low'); % use 4th order butterworth to further reduce hf noise
rawtrace = [0,diff(filter(b,a,rawtrace))];
thresh = 2*std(rawtrace); % Use 2 stdev threshold

binlicks = zeros(1,length(frames));

n=200; % Discard 200 datapoints (ringing)
while n < length(rawtrace) 
    if rawtrace(n)>thresh
        binlicks(n)=1;
        n=n+49; % ca. 50 ms interval for 1kHz data
    else
        n=n+1;
    end
end

for f = max(frames):-1:1
    lickcount(f) = sum(binlicks(frames == f));
end

end