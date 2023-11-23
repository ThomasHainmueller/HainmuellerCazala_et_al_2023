function trace = ExtractRoi(array, roi)
% Extracts the mean content (trace) of a roi (2D-logical Mask) in a 3D array
[x,y]=find(roi);

for n = 1:length(x)
    traces(n,1:size(array,3))=array(x(n),y(n),:);
end
trace = mean(traces,1);
end