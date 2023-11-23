function trace = remove_50Hz(trace,fs)
% Remove 50 Hz and harmonics from a given trace.

fn=fs/2; % Nyquist frequency
f0 = 50:50:fn; % For 50 Hz and harmonics
freqRatio = f0/fn; 
notchWidth = 0.1;

for n=1:length(f0)
    notchZeros = [exp( sqrt(-1)*pi*freqRatio(n) ), exp( -sqrt(-1)*pi*freqRatio(n) )];
    notchPoles = (1-notchWidth) * notchZeros;
    b = poly( notchZeros ); %# Get moving average filter coefficients
    a = poly( notchPoles ); %# Get autoregressive filter coefficients
    trace = filter(b,a,trace);
end;

end