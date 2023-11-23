function [rois]=CCRois(centers,array)
% extracts a set of rois around the points given in 'centers'. The rois
% include all points within a certain radius whose signals are sufficiently
% correlated to the central signal.
% centers: logical, ROI centers 
% array: 3D double, imaging video

% Sigma and dimensions for the gaussian kernel of the ROI determinate.
sigma = 4;
KerDim = 10;

% Threshold for the weighted cross correlation
CCthreshold = 0.1;

% Generates a matrix 'kernel' with a 2D gaussian probability distribution 
[xgrid, ygrid] = meshgrid(-KerDim:KerDim);
kernel = normpdf(xgrid,0,sigma).*normpdf(ygrid,0,sigma);
kernel = kernel./kernel(KerDim+1,KerDim+1);

% Find indices of centers in array for linear indexing
%indices = find(centers);
% Find indices of centers
[x, y] = find(centers);

% expand array so that rois with centers close to border don't 
% exceed array bounds during calculation
array = padarray(array,[(KerDim+1) (KerDim+1) 0],NaN);
rois = zeros(size(array,1),size(array,2),length(x));
x = x+KerDim+1;
y = y+KerDim+1;

for n = length(x):-1:1
    % region = array(x(n)-(size+1):x(n)+size,y(n)-(size+1):y(n)+size);
    % localCC = zeros(size(kernel,1),size(kernel,2));
    % Center pixel
    % thing1 = reshape(array(x(n),y(n),:)-mean(tc(y,x,:),3),[1 1 numFrames]); % Extract center pixel's time course and subtract its mean
    center = single(reshape(array(x(n),y(n),:)-mean(array(x(n),y(n),:)),[1 1 size(array,3)]));
    ad_a   = single(sum(center.*center,3));    % Auto corr, for normalization later
        
    % Neighborhood
    a = single(array(x(n)-KerDim:x(n)+KerDim,y(n)-KerDim:y(n)+KerDim,:));       % Extract the neighborhood
    b = single(mean(array(x(n)-KerDim:x(n)+KerDim,y(n)-KerDim:y(n)+KerDim,:)));  % Get its mean
    surround = bsxfun(@minus,a,b);                  % Subtract its mean
    ad_b = single(sum(surround.*surround,3));               % Auto corr, for normalization later
        
    % Cross corr
    localCC = sum(bsxfun(@times,center,surround),3)./sqrt(bsxfun(@times,ad_a,ad_b)); % Cross corr with normalization
    weighted = localCC.*kernel;
    ROI = weighted > CCthreshold;
    rois(x(n)-KerDim:x(n)+KerDim,y(n)-KerDim:y(n)+KerDim,n)=ROI;
end
% Remove the parts that were padded for border correction
rois = rois(KerDim+2:end-(KerDim+1),KerDim+2:end-(KerDim+1),:);
rois = logical(rois);
end