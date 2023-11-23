function [category] = numcategory(value,divisions,maxvalue)
% Classifies a value as beein a member of the nth division of a value range
% extending to maxvalue. E.g. divisions = 100 would yield a perecentage
% rank of the value.

if value > maxvalue
    category = divisions;
end

maxvalue = maxvalue/divisions;

for n = divisions:-1:1
    if value<=n*maxvalue && value>(n-1)*maxvalue
        category = n;
    end
end

end