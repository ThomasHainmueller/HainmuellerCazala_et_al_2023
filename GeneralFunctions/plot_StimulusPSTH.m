function plot_StimulusPSTH(data,cellno,varargin)
% Plot the activity of one cell in a dataset of the new type for all
% categories in detail (i.e. for every run). Good to check trace detection
% and place-cell identification.

args=struct('category',3,'PSTH_window',-90:150,'subst_baseline',true,'stim_onset',90);

% Overwrite default parameters if required
for pair = reshape(varargin,2,[])
    if isfield(args,pair{1})
        args.(pair{1})=pair{2};
    else
        error('Input argument does not exist\n')
    end
end


figure
cols = 5; % four stimuli: Visual, Audio, airpuff, reward, speed
titles = {'visual','audio','airpuff','reward'};

% Iterate through the stimuli
for s = 1:4
    
    PSTH = data.cells{cellno}.categories{args.category}.StimulusPSTH{s};
    
    if args.subst_baseline
        for r=1:size(PSTH,1)
            PSTH(r,:)=PSTH(r,:)-nanmean(PSTH(r,1:args.stim_onset));
        end
    end
    
    subplot(2,cols,s)
    % Mean +/- std for dF over y coord for all runs together
    hold on
    plot(args.PSTH_window,nanmean(PSTH),'k')
    plot(args.PSTH_window,nanmean(PSTH)-nanstd(PSTH),'r')
    plot(args.PSTH_window,nanmean(PSTH)+nanstd(PSTH),'r')
    title(titles{s});
    ylim([-.5 1]);
    hold off

    subplot(2,cols,cols+s)
    % dF over y coord for each run in this category
    imagesc(PSTH);
    colormap('jet')
    % TODO: reset to -0.5 1.0
    caxis([-0.5 1.0])
end

% Get dF over speed
thisdFodY = deal(data.cells{cellno}.categories{args.category}.dFodY);
thisdFodY = cat(2,thisdFodY{:});
thisdFodY = thisdFodY';

subplot(2,cols,5)
hold on
% Mean +/- std for dF over speed for all runs together
plot(nanmean(thisdFodY),'k')
plot(nanmean(thisdFodY)-nanstd(thisdFodY),'r')
plot(nanmean(thisdFodY)+nanstd(thisdFodY),'r')
title('speed');
hold off

subplot(2,cols,10)
% dF over speed
imagesc(thisdFodY);
colormap('jet')
caxis([-0.5 1.0])

end
%     subplot(2,cols,[1 cols+1])   
%     % All traces on top of each other, each as y coord / trace pair
%     if size(thisdFoT,1) > size(thisdFoT,2)
%         thisdFoT = thisdFoT'; % 171209 correction
%         thisft = thisft';
%         thisy = thisy';
%         if significantonly
%             thismask = thismask';
%         end
%     end
%         for iter = size(thisdFoT,1):-1:1
%             hold on;
%             plot(thisft(iter,:),thisdFoT(iter,:)-2*iter,'Color',[.5 .5 .5]);
%             if significantonly
%                 plot(thisft(iter,:),thismask(iter,:)-2*iter,'r');
%             end
%             plot(thisft(iter,:),thisy(iter,:)-2*iter+1);
%         end
%         title(['Category ' num2str(c)])
% 
%     if dFodY
%         subplot(2,cols,3)
%             % dF over speed for each run in this category
%             imagesc(thisdFodY);
%             colormap('jet')
%             caxis([-0.1 1.0])
%         subplot(2,cols,6)
%             % Mean +/- std for dF over y coord for all runs together
%             plot(nanmean(thisdFodY),'k')
%             hold on
%             plot(nanmean(thisdFodY)-nanstd(thisdFodY),'r')
%             hold on
%             plot(nanmean(thisdFodY)+nanstd(thisdFodY),'r')
%     end
% end
% end

