function [normtrace, background] = normtobaseline(trace, Fbaseline)
% Normalize a trace to the 'Fbaseline' conjunctive % of it's lowest values.


if nargin < 2
    Fbaseline = 0.1;
end

window = round(length(trace)-Fbaseline*length(trace));

for n = length(trace)-window:-1:1
    background(n) = mean(trace(n:n+window));
end

background = min(background);

normtrace = (trace-background)/background;

clear trace
end