function bounds=bidirectional_bounds(image)
% Find the first and last pixel (dim 1) of a bidirectional scan image where
% the lines are not overlapping (e.g. the usefull area) when put in the 2nd
% dimension of the image (e.g. image(:,first:last))

thresholdfactor = 0.1; % i.e. minimum + 10% is considered background brightness

forwardhist = mean(image(1:2:end,:),1);
reversehist = mean(image(2:2:end,:),1);

if any(reversehist>65500)
    rf = find(reversehist>65500,1,'last');
    reversehist(reversehist>65500)=min(reversehist);
end

forwardthresh = min(forwardhist) + thresholdfactor*(max(forwardhist)-min(forwardhist));
reversethresh = min(reversehist) + thresholdfactor*(max(reversehist)-min(reversehist));

ff = find(forwardhist>forwardthresh,1,'first');
fl = find(forwardhist>forwardthresh,1,'last');

if ~exist('rf','var')
    rf = find(reversehist>reversethresh,1,'first');
end

rl = find(reversehist>reversethresh,1,'last');

bounds(1) = max([ff,rf]);
bounds(2) = min([fl,rl]);
end
% threshold = 1E7; % Delta intensity for vertical that is considered unusually large 
% band = 10; % Safety margin (in pixels) that is certainly not part of the image
% 
% % Loop through vertical image lines and calculate 1st differential
% for n=size(image,2):-1:1
%     t(:,n)=diff(image(:,n)).^2;
% end
% 
% t = mean(t,1);
% t = t(band:end);
% 
% first = find(t<threshold,1,'first')+band;
% if find(t(first:end)>threshold,1,'first')
%     last = find(t(first:end)>threshold,1,'first')+first;
% else
%     last = size(t,2); % Set to end of image
% end
% %last=find(t(first:end)<threshold,1,'last')+first;
% end