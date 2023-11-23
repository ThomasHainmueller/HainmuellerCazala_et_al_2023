function data = salvage_stimuli(data, category, varargin)
% Function to correct the bug in the acquisition of x- and y-traces
% associated with the airpuffs present in the datasets up until 180927.

args=struct('transients',false,'stim_amplitudes',[.5, 1, 1.5, 2],...
    'stim_range',.1,'stim_mindur',10,'stim_length',30,'stim_interval',126); 
    % transients: Use significant transients only (n.A. for interneurons)
    % stim_amplitudes: Amplitude of the x-transients set in VR
    % stim_range: Amplidude window around each stim level
    % stim_mindur: Minimum duration, to avoid missclassification of slopes
    % stim_length, stim_interval: In frames, for the reconstructed y trace.
    % 15 and 63 for unidirectional, 30 and 126 for bidirectional

% Overwrite default parameters if required
for pair = reshape(varargin,2,[])
    %args.(pair{1}) = pair{2};
    if isfield(args,pair{1})
        args.(pair{1})=pair{2};
    else
        error('Input argument does not exist\n')
    end
end

for r = 1:length(data.metadata.categories{category}.x)
    thisx = data.metadata.categories{category}.x{r};
    thisy = data.metadata.categories{category}.y{r};
    tracelength = length(thisx);
    
    thisx = thisx(~isnan(thisx));
    thisy = thisy(~isnan(thisy));
    
    stimuli = classify_x(thisx, args.stim_amplitudes, args.stim_range, args.stim_mindur);
    
    airpufftimes = find([0,diff(stimuli(3,:))==1]); % Assuming airpuffs are the 3rd stimulus
    
    % get onsettimes of the remaining stimuli.
    onsettimes = [];
    for s = [1,2,4]
        onsettimes = [onsettimes, find([0,diff(stimuli(s,:))==1])];
    end
    onsettimes = sort(onsettimes);
    
    keepx = false(1,length(thisx));
    %keepx(1:airpufftimes(1)-1) = true;
    %keepx = false(1,tracelength);
    keepx(1:airpufftimes+args.stim_length+args.stim_interval)=true;
    
    for ap = 2:length(airpufftimes)
        if find(onsettimes>airpufftimes(ap-1),1)
            start = onsettimes(find(onsettimes>airpufftimes(ap-1),1));
            %keepx(start:airpufftimes(ap)-1) = true;
            keepx(start:airpufftimes(ap)+args.stim_length+args.stim_interval)=true;
        else
            break
        end
    end
    
    if find(onsettimes>airpufftimes(end),1)
        start = onsettimes(find(onsettimes>airpufftimes(end),1));
        keepx(start:end) = 1;
    end
    
    if length(keepx)>length(thisx)
        keepx = keepx(1:length(thisx));
        %keepx = keepx(1:tracelength);
    end
    figure; hold on; plot(thisx); plot(keepx);
    thisx = thisx(keepx);
    thisx(end+1:tracelength) = 0;
    thisy = thisy(keepx);
    thisy(end+1:tracelength) = 0;
    
    data.metadata.categories{category}.x{r} = thisx;
    data.metadata.categories{category}.y{r} = thisy;
    
%     % fill new x trace
%     newx(1:airpufftimes(1)-1) = thisx(1:airpufftimes(1)-1);
%     for ap = 1:length(airpufftimes)-1
%         newx(end+1:end+args.stim_length)=1.5; % Assuming airpuff voltage at 1.5V
%         newx(end+1:end+args.stim_interval)=0;
%         nextstim = onsettimes(find(onsettimes>length(newx),1));
%         newx = [newx,thisx(nextstim:airpufftimes(ap+1)-1)];
%     end
%     newx(end+1:end+args.stim_length)=1.5;
%     newx(end+1:end+args.stim_interval)=0;
%     
%     if length(newx)>length(thisx)
%         newx = newx(1:length(thisx));
%     else
%         nextstim = onsettimes(find(onsettimes>length(newx),1));
%         newx = [newx,thisx(nextstim:length(thisx))];
%         newx(end+1:length(thisx)) = 0;
%     end
        
end
end

% % get stimulus timestamps
% xdata = cat(2, data.metadata.categories{category}.x{:});
% stimuli = classify_x(xdata, args.stim_amplitudes, args.stim_range, args.stim_mindur);
% 
% for n=1:length(data.cells)
%     thistrace = cat(2, data.cells{n}.categories{category}.dFoT{:});
%     for s=1:size(stimuli,1)
%         onsettimes = find([0,diff(stimuli(s,:))==1]);
% 
%         for t = 1:length(onsettimes)
%             if onsettimes(t)+args.PSTH_wdw(2)<length(thistrace) && onsettimes(t)-args.PSTH_wdw(1)>0
%                 thisPSTH(t,:) = thistrace(onsettimes(t)-args.PSTH_wdw(1):...
%                     onsettimes(t)+args.PSTH_wdw(2));
%             end
%         end
%         %meanPSTH(s,:) = mean(thisPSTH,1); Modified 180920
%         PSTH{s} = thisPSTH;
%         clear thisPSTH; % CLEAR variable to avoid inheritance between stimuli
%     end
%     %data.cells{n}.categories{category}.StimulusPSTH = meanPSTH; 180920
%     data.cells{n}.categories{category}.StimulusPSTH = PSTH;
% end
% end

function stimuli = classify_x(xdata, stim_amplitudes, stim_range, stim_mindur)
% Inner function, returns a logical nstimuli x nframes matrix
% indicating when each of the stimuli is active.
for s = 1:length(stim_amplitudes)
    thisstim = and(xdata>stim_amplitudes(s)-stim_range, xdata<stim_amplitudes(s)+stim_range);
    thisstim = medfilt1(single(thisstim), stim_mindur);
    stimuli(s,:) = logical(thisstim);
end
end