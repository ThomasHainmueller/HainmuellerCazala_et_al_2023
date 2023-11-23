function stimplot = stimulus_psth_meta(dataset,category,varargin)
% Plots lineartrack acquired activity grouped to category for multiple
% datasets. Takes dFoY traces as given in the dataset - should allready be
% recalculated only for the significant transients. Parameters (activity,
% SI or PF p-values) are checked for the given categories and cells
% meeting the criteria will be plotted in an 'and' or 'or' manner as
% specified. Generate appropriate input by placing several tdata
% (derived from get_lineartracksX => convert_data) into a cell array.

if nargin < 2
    category = 3;
end

args=struct('signal','zscored','sortstimulus',1,'plotmean',true,'response_wdw',60,...
    'stim_onset',90,'plotrange',[-.1,.4]); 
    % SI/PF significance levels of spatial info/placefields from bootstrap
    % sortstimulus: 1=visual, 2=audio, 3=airpuff, 4=reward
    % response window: for stimlus response estimation, in # frames
    % stim onset: first frame (from beginning of the stored PSTH) at which
    %   the stimuli begin.

% Overwrite default parameters if required
for pair = reshape(varargin,2,[])
    %args.(pair{1}) = pair{2};
    if isfield(args,pair{1})
        args.(pair{1})=pair{2};
    else
        error('Input argument does not exist\n')
    end
end

labels = {'visual','audio','airpuff','reward'};
stimplot = [];
respindex = []; % response indices

% Get Stimulus PSTH for all cells throughout datasets.
for d=1:numel(dataset)
    ds=dataset{d};
    
    for n=1:length(ds.cells)
        for s=1:length(ds.cells{n}.categories{category}.StimulusPSTH)
            thisPSTH = nanmean(ds.cells{n}.categories{category}.StimulusPSTH{s});
            %thisstimplot(s,:) = thisPSTH;
            thisstimplot(s,:) = thisPSTH-nanmean(thisPSTH(1:args.stim_onset)); % Adjust baseline levels
            thisrespindex(s) = nanmean(...
                thisPSTH(args.stim_onset:args.stim_onset+args.response_wdw))/...
                nanmean(thisPSTH(args.stim_onset-args.response_wdw:args.stim_onset-1));
        end
        respindex(end+1,:) = thisrespindex;
        stimplot(:,:,end+1) = thisstimplot; % increment only once per cell!
        clear thisstimplot;
    end
end

stimplot(:,:,1) = []; % Remove the first (empty) slice

[~,sortindices] = sort(respindex(:,args.sortstimulus),'descend');

% figure; plot(respindex(:,args.sortstimulus));
% figure; imagesc(squeeze(stimplot(args.sortstimulus,:,sortindices))'); title('sorted')
% figure; plot(squeeze(nanmean(stimplot(args.sortstimulus,:,:),3)));
% return

% Display part
figure;
for s = 1:size(stimplot,1)
    subplot(2,size(stimplot,1),s)
        imagesc(squeeze(stimplot(s,:,sortindices))');
        title(labels{s});
        colormap(jet)
        caxis(args.plotrange);
    subplot(2,size(stimplot,1),s+size(stimplot,1))
        hold on
        plot(squeeze(nanmean(stimplot(s,:,:),3)),'k');
        plot(squeeze(nanmean(stimplot(s,:,:),3))-...
            (std(squeeze(stimplot(s,:,:)),0,2)/sqrt(size(stimplot,3)))','r');
        plot(squeeze(nanmean(stimplot(s,:,:),3))+...
            (std(squeeze(stimplot(s,:,:)),0,2)/sqrt(size(stimplot,3)))','r');
        ylim([-.1,.2])
end
end
