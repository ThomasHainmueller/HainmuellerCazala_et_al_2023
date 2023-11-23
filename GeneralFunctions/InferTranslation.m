function tform = InferTranslation(registered)
% Infers translation as the minimum number of full pixels shifted.
% Based on 'all-zero' rows / columns.

midpoint=[round(size(registered,1)/2),round(size(registered,2)/2)];
base_points=[midpoint;midpoint+1];
input_points=base_points;

ZeroRows = find(all(registered==0,2));
ZeroColumns = find(all(registered==0,1));

UpperZeros = length(find(ZeroRows<=midpoint(1)));
LowerZeros = length(find(ZeroRows>=midpoint(1)));
LeftZeros = length(find(ZeroColumns<=midpoint(2)));
RightZeros = length(find(ZeroColumns>=midpoint(1)));

if UpperZeros
    input_points(:,2)=input_points(:,2)-UpperZeros;
elseif LowerZeros
    input_points(:,2)=input_points(:,2)+LowerZeros;
end

if LeftZeros
    input_points(:,1)=input_points(:,1)-LeftZeros;
elseif RightZeros
    input_points(:,1)=input_points(:,1)+RightZeros;
end

tform = cp2tform(input_points,base_points,'nonreflective similarity');
end