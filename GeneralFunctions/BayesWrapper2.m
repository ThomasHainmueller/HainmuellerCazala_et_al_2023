function [pPost,pPostMax,DecErr,Y] = BayesWrapper2(data, varargin)
% Get traces from the selected categories of a dataset and perform bayesian
% decoding of the trajectory using the traces from specified cells
% Optional arguments:
% categories:	int, default: use all categories of the dataset
% cells:        int, specifiy the cells used for en- and decoding
% nsegments:	int, decoding is performed on segments omited during
%           	encoding. Specify the number of segements.
% moving_only:  bool, default true.
% transients_only:	bool, whether to use only detected transients or the
%               	entire fluorescence trace.
% spbins:       double, spatial bins for discretizing the spatial maps
% TauDecode:    double, size of the time bins for decoding in seconds; 
%               Default sampling frequency. Downsamples using resample()
% Track_length: double, length of the linear track, default 4 m, to convert
%               decoding error in m.
% savefigure:    bool, whether to save the figure after displaying.
% maxNruns:     int, maximum number of runs used to decode. Will either
%               this number or the maximum available in the smallest
%               category. Set to 'Inf' to decode all runs.
%   
% Output:
% pPost:        {nruns x ncategories} cell of [ntimebins x nspatialbins x
%               ncategories] matrices denoting the decoded probability (or
%               correlation value, for the PV decoder).
% pPostMax:     {nruns x ncategories} cell of [ntimebins x 2]
%               matrices holding the most likely decoded position (:,1) and
%               category (:,2) for each timestamp.
% DecErr:       {nruns x ncategories} cell of [ntimebins x 2]
%               matrices holding the position error in meters (:,1) and
%               category error (:,2) for each timebin relative to true
%               location. NOTE: If more than two categories are compared,
%               this should be adjusted unless they have an interval
%               relationship.
%
% Note to self: Positional- and Contextual information are the respective
% marginals of the two decoded probability distributions.
% 
% TODO: 
% - Omit cells with high zero-timelag correlations (same cell twice) on the
%   dFoT traces (for the gradual exclusion at least).
% - Simple way to get chance-levels of decoding: shuffle cell IDs in the
%   template after encoding.
%
%% Input handling
p = inputParser;
addParameter(p,'categories',[],@isnumeric);
addParameter(p,'cells',[],@isnumeric);
addParameter(p,'nsegments',10,@isnumeric);
addParameter(p,'moving_only',true,@islogical);
addParameter(p,'transients_only',false,@islogical);
addParameter(p,'spbins',0.1:0.025:2.1,@isnumeric);
addParameter(p,'TauDecode',0,@isnumeric);
addParameter(p,'display',false,@islogical);
addParameter(p,'savefigure',false,@islogical);
addParameter(p,'Track_length',4,@isnumeric);
addParameter(p,'maxNruns',20,@isnumeric);

parse(p,varargin{:})

nsegments = p.Results.nsegments;
moving_only = p.Results.moving_only;
transients_only = p.Results.transients_only;
spbins = p.Results.spbins;
display = p.Results.display;
savefigure = p.Results.savefigure;
Track_length = p.Results.Track_length;
maxNruns = p.Results.maxNruns;

% Compatibility issue
if iscell(data.metadata.categories{1}.acquisition_rate)
    acquisition_rate = data.metadata.categories{1}.acquisition_rate{1};
else
    acquisition_rate = data.metadata.categories{1}.acquisition_rate(1);
end

if ~isempty(p.Results.categories)
    categories = p.Results.categories;
else
    categories = 1:length(data.metadata.categories);
end

if ~isempty(p.Results.cells)
    cells = p.Results.cells;
else
    cells = 1:length(data.cells);
end

if p.Results.TauDecode ~= 0
    TauDecode = p.Results.TauDecode;
else
    TauDecode = 1/acquisition_rate;
end

% Get number of runs for decoding (should be same for all categories)
for c = categories
    nruns(c) = length(data.metadata.categories{c}.y);
end

nruns = min([nruns(nruns>0) maxNruns]);

%% Iterate through runs, get dFoY for each run for building templates
for r = nruns:-1:1
    for c = flip(categories)
        
        % Get the mouse's position on the track
        if moving_only
            thisy = data.metadata.categories{c}.y{r}.*...
                data.metadata.categories{c}.moving{r};
        else
            thisy = data.metadata.categories{c}.y{r};
        end
        
        % Get the calcium traces
        for n = length(data.cells):-1:1
            if transients_only
                thistr = data.cells{n}.categories{c}.dFoT{r}.*...
                    data.cells{n}.categories{c}.transientmask{r};
            else
                thistr = data.cells{n}.categories{c}.dFoT{r};
            end
            
            % Trace should be all >= 0. TRY DIFFERENT FIXES HERE
            thistr(thistr<0) = 0;
            
            dFoYs(:,r,c,n) = SBdiscretize(thistr,thisy,spbins);
            % MEMO: Remove NaN values after creating maps!
        end
    end
end

if display
    h = figure('Position',[10 50 1900 900]);
    try
        SessName = split(data.metadata.categories{1}.filename{1},'_');
    catch
        SessName{1} = 'DateUnknown'; SessName{2} = 'AnimalUnknown';
    end
    sgtitle(sprintf('%s - %s - %d cells %.2f seconds',SessName{1},SessName{2},length(cells),TauDecode));
end

%% Encoding and decoding part
for r = nruns:-1:1
    % Create template
    template_runs = true(1,nruns);
    template_runs(r) = false;
    
    template = nanmean(dFoYs(:,template_runs,:,:),2);
    template = reshape(template,size(template,1)*size(template,3),size(template,4));
    template(isnan(template))=0;
    
    for c = flip(categories)
        % Retrieve transients for the run to be decoded
        for n = length(data.cells):-1:1
            if transients_only
                thistr = data.cells{n}.categories{c}.dFoT{r}.*...
                    data.cells{n}.categories{c}.transientmask{r};
            else
                thistr = data.cells{n}.categories{c}.dFoT{r};
            end
            % Set NaN to 0; may be problematic for Baeysian approach
            thistr(isnan(thistr))=0;
            
            % No values < 0; TRY DIFFERENT FIXES HERE; SEE L. 72
            thistr(thistr<0)=0;
            thistr = double(thistr);
            
            % Downsampling
            %nTs = round(1/data.metadata.categories{c}.acquisition_rate(r)/TauDecode*length(thistr));
            nTs = round(1/acquisition_rate/TauDecode*length(thistr)); % Assume homogenous aquisition rate throughout.
            
            decode_tr(n,:) = resample(thistr,nTs,length(thistr));
        end
        
        % Decoding
        % pPost{r,c} = bayes_decode(template,thistr,1/data.metadata.categories{c}.acquisition_rate(r));
        [~, pPost{r,c}] = PVcorr_decode(template,decode_tr);
        [~,thispos] = max(pPost{r,c},[],2);
        
        if moving_only
            moving = data.metadata.categories{c}.moving{r};
            moving = resample(double(moving),nTs,length(moving));
            moving = moving>max(moving)/2; % re-binarize
            thispos(~moving) = NaN;
        end
        
        pPostMax{r,c}(:,1) = mod(thispos,size(pPost{r,c},2)/length(categories)); % Get position
        pPostMax{r,c}(:,2) = ceil(thispos/size(pPost{r,c},2)*length(categories)); % Get context
        clear decode_tr
        
        
        % Calculate decoding error
        thisy = (data.metadata.categories{c}.y{r}-spbins(1))./(spbins(end)-spbins(1))*(length(spbins)-1);
        thisy(thisy<0) = 0;
        thisy(thisy>(length(spbins)-1)) = length(spbins)-1;
        thisy = resample(thisy,nTs,length(thisy));          
        
        % Plotting results if desired
        if display
            %subplot(length(categories),nruns,r+(c-1)*nruns)
            %axpos = get(gca, 'Position');
            axpos(1) = (r-1)/nruns+.01; % x
            axpos(2) = 1.02-(c/length(categories)); %- 1/length(categories) +.02; %y
            axpos(3) = 1/nruns-.015; % width
            axpos(4) = 1/length(categories)-.1; % height
            %subplot(length(categories),nruns,r+(c-1)*nruns,'Position',axpos)
            axes('Position',axpos);
            %set(gca, 'Position', axpos)  
            hold on;
            
            ts = linspace(0,length(thisy)*TauDecode,length(thisy));
            
            imagesc([min(ts) max(ts)], [1 size(pPost{r,c},2)],pPost{r,c}'); 
            colormap(flipud(gray));
            plot(ts,thisy+(c-1)*(length(spbins)-1),'r');
            
            scatter(ts,thispos,3,'b');
            
            ylim([1 size(pPost{r,c},2)]);
            xlim([1 max(ts)]);
 
        end
        
        pPost{r,c} = reshape(pPost{r,c},size(pPost{r,c},1),...
            size(pPost{r,c},2)/length(categories),length(categories));
        
        DecErr{r,c}(:,1) = abs(pPostMax{r,c}(:,1)-thisy')/...
            size(pPost{r,c},2)*Track_length;
        DecErr{r,c}(:,2) = abs(pPostMax{r,c}(:,2)-c);
        
        Y{r,c} = thisy';
    end
end

if savefigure
    saveas(h,sprintf('%s-%s-%d_cells-%.2f_seconds.png',SessName{1},SessName{2},length(cells),TauDecode));
end

end
