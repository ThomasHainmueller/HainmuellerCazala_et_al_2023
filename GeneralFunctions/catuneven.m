function [a] = catuneven(elements, pad)
% Function to concatenate vectors of uneven size.
% Concatenates the elements given in the cellarray and pads
% empty spaces with the element specified in pad (e.g. zero, NaN).
% 
% New feature 11/21/2023 - Allows for concatenating multiple arrays, by
% default along second (i.e. column) dimension.
% 
% Help: Row vector = ones(1,10), 1 row, 10 columns
%       Column vector = ones(10,1), 10 rows, 1 column

le = []; nel = [];

for n = 1:length(elements)
    
    % Force column vectors if singular, for arrays assume that observations
    % are in rows, variables in columns
    if isrow(elements{n})
        elements{n} = elements{n}';
    end
    
    le = [le size(elements{n},1)];
    
    % Preserve structure of input if empty element given
    thisnel = max(size(elements{n},2),1);
    nel = [nel thisnel];
end

maxle = max(le);

%a = NaN(maxle,length(elements));
a = NaN(maxle,sum(nel));
a(:) = pad;

for n = 1:length(elements)
    v = elements{n};
    if n>1
        i1 = sum(nel(1:n-1))+1;
    else
        i1 = 1;
    end
    a(1:size(v,1),i1:i1+nel(n)-1)=v;
    %a(i1:i1+nel(n)-1,1:size(v,2));
end

% for n = length(elements):-1:1
%     v = elements{n};
%     %a(n,1:length(v))=v;
%     a(1:length(v),n)=v;
% end
%a=a';
end