function data = spatial_param_rnd(data,varargin)
% Calculate Moving/Immobile ratio, speed-modulation slope, SI, SI/AUC, 
% trial-trial variability, 1st to 2nd half stability for real- and
% circularly shuffled data.
p = inputParser;
addParameter(p,'nshuffles',100,@isnumeric);
addParameter(p,'bins',.1:.025:1.9,@isnumeric);
addParameter(p,'category',1,@isnumeric);

parse(p,varargin{:})

nshuffles = p.Results.nshuffles;
bins = p.Results.bins;
ca = p.Results.category;

rng(1);
framerate = data.metadata.categories{1}.acquisition_rate(1);

%% Generate LUT for extracting per-bin information
mvg = cat(2,data.metadata.categories{ca}.moving{:});
y = cat(2,data.metadata.categories{ca}.y{:});
y = y(mvg);

% Create logical matrix to discretize position
for b = length(bins):-1:2
    yb(b-1,:) = y>bins(b-1) & y<=bins(b);
end

ryb = false(size(yb,1),size(yb,2),nshuffles);
ry = zeros(length(y),nshuffles);
%for sh = nshuffles:-1:1
parfor sh = 1:nshuffles
    ry(:,sh) = circshift(y,round(rand*length(y)));
    %ryb(:,:,sh) = circshift(yb,round(rand*size(yb,2)),2);
end

%% Calculate random spatial maps and extract information
for n = 1:length(data.cells)
    % TODO: Add p-value and first shuffle value to dataset for display.
    tr = cat(2,data.cells{n}.categories{ca}.dFoT{:});
    tr = tr(mvg);
    
%     tr = repmat(tr,size(yb,1),1);
    
    mid = round(size(tr,2)/2);
    
%     map = tr.*yb;   
%     map(~yb) = NaN; % Remove zero values
%     dFoY = nanmean(map,2);
    dFoY = discretize(tr,y,bins);
    % Create separate mean traces for 1st and 2nd half
%     dFoY1 = nanmean(map(:,1:mid),2);
%     dFoY2 = nanmean(map(:,mid:end),2);
    dFoY1 = discretize(tr(1:mid),y(1:mid),bins);
    dFoY2 = discretize(tr(mid:end),y(mid:end),bins);
    SI = spatial_info(tr,ry(:,sh)',bins,framerate);
    
%     for sh = nshuffles:-1:1
%         for b = 1:size(ryb,1)
%             dFoYr(sh,b) = nanmean(tr(ryb(b,:,sh)));
%             dFoYr1(sh,b) = nanmean(tr(ryb(b,1:mid,sh)));
%             dFoYr2(sh,b) = nanmean(tr(ryb(b,mid:end,sh)));
%         end
%     end
    
%     dFoYr = zeros(nshuffles,size(dFoY,1));
%     dFoYr1 = zeros(nshuffles,size(dFoY,1));
%     dFoYr2 = zeros(nshuffles,size(dFoY,1));
%     for sh = nshuffles:-1:1
    for sh = 1:nshuffles
%         rmap = tr.*ryb(:,:,sh);
%         rmap(~ryb(:,:,sh)) = NaN;
%         dFoYr(sh,:) = nanmean(rmap,2);
%         dFoYr1(sh,:) = nanmean(rmap(:,1:mid),2);
%         dFoYr2(sh,:) = nanmean(rmap(:,mid:end),2);
        dFoYr(sh,:) = discretize(tr,ry(:,sh),bins);
        dFoYr1(sh,:) = discretize(tr(1:mid),ry(1:mid,sh),bins);
        dFoYr2(sh,:) = discretize(tr(mid:end),ry(mid:end,sh),bins);
        
       SIr(sh) = spatial_info(tr,ry(:,sh)',bins,framerate);
    end
%     figure; imagesc(dFoYr);
%     hold on; plot(dFoY*1000,'k');
%     plot(nanmean(dFoYr,1)*10000);
% 
%     figure; hold on;
%     plot(nanmean(dFoYr,1),'color',[.5 .5 .5]);
%     plot(dFoY,'k'); 
%     plot(dFoY1,'r');
%     plot(dFoY2,'b');
end
end