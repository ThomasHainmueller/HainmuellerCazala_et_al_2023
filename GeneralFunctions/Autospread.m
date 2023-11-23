% Automatism to get the GRspread vector of an subindex in a 'mestaghandle' file

function peakamplitudes = Autospread(file, index, fFilter, SpAvg, F0start, F0end, tRecord)

if nargin < 7
    tRecord = 4;
end 

if nargin < 6
    F0end = 1.2;
end

if nargin < 5
    F0start = 0.3;
end

if nargin < 4
    SpAvg = 1;
end

if nargin < 3
    fFilter = 100;
end

% main executive part of the function
x = getGR(file,[index]);
% extract matrix from cell
x = x{1};
% x is a GoR matrix of type double (x,y) now
x = GoR(x,1,1,F0start, F0end, tRecord);
peakamplitudes = GRspread(x,F0start,tRecord,fFilter,SpAvg);
end