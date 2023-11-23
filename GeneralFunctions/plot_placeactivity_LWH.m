function plot_placeactivity_LWH(data,cellno,significantonly,dFodY)
% Plot the activity of one cell in a dataset of the new type for all
% categories in detail (i.e. for every run). Good to check trace detection
% and place-cell identification.
if nargin < 4
    dFodY = false;
end
if nargin < 3
    significantonly = false;
end

for c = 1:length(data.cells{cellno}.categories)
%for c = 1;
    % Evertything for the individual trace plotting
    thisdFoT = deal(data.cells{cellno}.categories{c}.dFoT);
    thisdFoT = catuneven(thisdFoT,NaN);
    thisy = deal(data.metadata.categories{c}.y);
    thisy = catuneven(thisy,NaN);
    thisft = deal(data.metadata.categories{c}.ft);
    thisft = catuneven(thisft,NaN);
    if significantonly
        thismask = deal(data.cells{cellno}.categories{c}.transientmask);
        thismask = catuneven(thismask,false);
        thismask = thisdFoT.*thismask;
        thismask(thismask==0)=NaN;
    end
    
    if significantonly
        thisdFoY = deal(data.cells{cellno}.categories{c}.dFoY);
        thisdFoY = cat(2,thisdFoY{:});
    else
        thisdFoY=[];
        for catrun = 1:length(data.cells{cellno}.categories{c}.dFoY)
            thisdFoY = cat(2,discretize(...
                data.cells{cellno}.categories{c}.dFoT{catrun},...
                data.metadata.categories{c}.y{catrun},0.1:0.05:2.1),thisdFoY);
        end
    end
    thisdFoY = thisdFoY';
    %thisdFoY(thisdFoY==0) = NaN;
    if dFodY
        thisdFodY = deal(data.cells{cellno}.categories{c}.dFodY);
        thisdFodY = cat(2,thisdFodY{:});
        thisdFodY = thisdFodY';
    end
    
    % ---- PLOTTING PART ----
    figure
    if dFodY
        cols=3;
    else
        cols=2;
    end
    subplot(2,cols,[1 cols+1])   
    % All traces on top of each other, each as y coord / trace pair
    if size(thisdFoT,1) > size(thisdFoT,2)
        thisdFoT = thisdFoT'; % 171209 correction
        thisft = thisft';
        thisy = thisy';
        if significantonly
            thismask = thismask';
        end
    end
        for iter = size(thisdFoT,1):-1:1
            hold on;
            plot(thisft(iter,:),thisdFoT(iter,:)-2*iter,'Color',[.5 .5 .5]);
            if significantonly
                plot(thisft(iter,:),thismask(iter,:)-2*iter,'r');
            end
            plot(thisft(iter,:),thisy(iter,:)-2*iter+1);
        end
        title(['Category ' num2str(c)])
    subplot(2,cols,2)
    % dF over y coord for each run in this category
        imagesc(thisdFoY);
        colormap('jet')
        % TODO: reset to -0.5 1.0
        caxis([-0.1 1.0])
    subplot(2,cols,cols+2)
    % Mean +/- std for dF over y coord for all runs together
        plot(nanmean(thisdFoY),'k')
        hold on
        plot(nanmean(thisdFoY)-nanstd(thisdFoY),'r')
        hold on
        plot(nanmean(thisdFoY)+nanstd(thisdFoY),'r')
    if dFodY
        subplot(2,cols,3)
            % dF over speed for each run in this category
            imagesc(thisdFodY);
            colormap('jet')
            caxis([-0.1 1.0])
        subplot(2,cols,6)
            % Mean +/- std for dF over y coord for all runs together
            plot(nanmean(thisdFodY),'k')
            hold on
            plot(nanmean(thisdFodY)-nanstd(thisdFodY),'r')
            hold on
            plot(nanmean(thisdFodY)+nanstd(thisdFodY),'r')

end
end

