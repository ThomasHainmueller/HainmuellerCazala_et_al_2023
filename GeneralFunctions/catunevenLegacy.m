function [a] = catunevenLegacy(elements, pad)
% Function to concatenate vectors of uneven size.
% Concatenates the elements given in the cellarray and pads
% empty spaces with the element specified in pad (e.g. zero, NaN).
%
% Help: Row vector = ones(1,10), 1 row, 10 columns
%       Column vector = ones(10,1), 10 rows, 1 column

le = [];

for n = length(elements):-1:1
    
    % Force column vectors, newly introduced 11/22/2023
    if isrow(elements{n})
        elements{n} = elements{n}';
    end
    
    le = [le size(elements{n},1)];
end

maxle = max(le);

%a = NaN(length(elements),maxle);
a = NaN(maxle,length(elements));
a(:) = pad;

for n = length(elements):-1:1
    v = elements{n};
    %a(n,1:length(v))=v;
    a(1:length(v),n)=v;
end
%a=a';
end