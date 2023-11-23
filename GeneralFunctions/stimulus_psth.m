function data = stimulus_psth(data, category, varargin)
% Takes a dataset and creates a PSTH for all recordings in one category for
% each cell in them. This assumes that different stimuli are labelled as
% transients with a distinct magnitude in the x-coordinate of the dataset.
% Useful for the IN-stimuli type of experiment to quantify responses to
% different sensory modalities.

args=struct('transients',false,'stim_amplitudes',[.5, 1, 1.5, 2],...
    'stim_range',.1,'stim_mindur',10, 'PSTH_wdw',[90,150]); 
    % transients: Use significant transients only (n.A. for interneurons)
    % stim_amplitudes: Amplitude of the x-transients set in VR
    % stim_range: Amplidude window around each stim level
    % stim_mindur: Minimum duration, to avoid missclassification of slopes
    % PSTH_wdw: Window (frames) around stimulus onset for PSTH calculation

% Overwrite default parameters if required
for pair = reshape(varargin,2,[])
    %args.(pair{1}) = pair{2};
    if isfield(args,pair{1})
        args.(pair{1})=pair{2};
    else
        error('Input argument does not exist\n')
    end
end

% get stimulus timestamps
xdata = cat(2, data.metadata.categories{category}.x{:});
stimuli = classify_x(xdata, args.stim_amplitudes, args.stim_range, args.stim_mindur);

for n=1:length(data.cells)
    thistrace = cat(2, data.cells{n}.categories{category}.dFoT{:});
    for s=1:size(stimuli,1)
        onsettimes = find([0,diff(stimuli(s,:))==1]);

        for t = 1:length(onsettimes)
            if onsettimes(t)+args.PSTH_wdw(2)<length(thistrace) && onsettimes(t)-args.PSTH_wdw(1)>0
                thisPSTH(t,:) = thistrace(onsettimes(t)-args.PSTH_wdw(1):...
                    onsettimes(t)+args.PSTH_wdw(2));
            end
        end
        %meanPSTH(s,:) = mean(thisPSTH,1); Modified 180920
        PSTH{s} = thisPSTH;
        clear thisPSTH; % CLEAR variable to avoid inheritance between stimuli
    end
    %data.cells{n}.categories{category}.StimulusPSTH = meanPSTH; 180920
    data.cells{n}.categories{category}.StimulusPSTH = PSTH;
end
end

function stimuli = classify_x(xdata, stim_amplitudes, stim_range, stim_mindur)
% Inner function, returns a logical nstimuli x nframes matrix
% indicating when each of the stimuli is active.
for s = 1:length(stim_amplitudes)
    thisstim = and(xdata>stim_amplitudes(s)-stim_range, xdata<stim_amplitudes(s)+stim_range);
    thisstim = medfilt1(single(thisstim), stim_mindur);
    stimuli(s,:) = logical(thisstim);
end
end