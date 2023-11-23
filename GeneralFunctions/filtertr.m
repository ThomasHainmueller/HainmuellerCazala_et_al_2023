% function to lowpass filter a cell of traces and return the peak of the filtered
% trace within a defined window. Uses Chebychev2 order 3 filtering
% fstop is the frequency, by which the attenuation first reaches 60dB, 20Hz
% have prooven well empirically, change if results don't look well
% wstart, wstop = window for peak detection (in seconds)
% fsamp = sampling frequency (in pts/second)

function [resfilt trfilt] = filtertr(trace, fstop, wstart, wstop, fsamp)

if nargin < 5
    fsamp = 100;
end

if nargin < 4
    wstop = 2.5;
end 

if nargin < 3 
    wstart = 1.2;
end

if nargin < 2
    fstop = 20;
end

fstopnorm = fstop/(fsamp/2);
wstart = wstart*fsamp;
wstop = wstop*fsamp;

% creating the digital filter
[z,p,k] = cheby2(3,60,fstopnorm);
[sos,g]=zp2sos(z,p,k);
Hd = dfilt.df2tsos(sos,g);

for c = [1:length(trace)]
    tr = trace(c).data;
    tr = reshape(tr,length(tr),1);
    trfilt(:,c) = filter(Hd,tr);
    resfilt(c) = max(trfilt(wstart:wstop,c));
    plot(trfilt(:,c));
    hold on;
end
resfilt = transpose(resfilt);
end
