function [segments] = FindEvents(metadata)
% Identify events (baseline, runs, memory) based on x,y traces and return
% a struct containing type, direction, start and end frames of the events.
% This function works with match-to-place behaviours recorded in the sequence 
% baseline - Sample run - delay period - test run - post behaviour.

frametimes = metadata{1};
xcord = metadata{2};
ycord = metadata{3};
dxcord=sqrt(diff(xcord).^2);
dycord=sqrt(diff(ycord).^2);

%figure,plot(xcord),hold on,plot(ycord),hold on,plot(dxcord),hold on,plot(dycord);

% X- or Y velocity threshold at which a trace start / end is assumed.
threshold = 0.4;

% Threshold for movement assumption
moving = 0.05;

% Indices at which transitions between segments might occur.
transitions = find(dxcord>threshold | dycord>threshold);
running_speed = dycord;
for n = [1:length(transitions)]
    running_speed(transitions(n)-1:transitions(n)+1)=0;
end

segments.baseline = [NaN,NaN];
segments.sample = [NaN,NaN];
segments.delay = [NaN,NaN];
segments.match = [NaN,NaN];
segments.post = [NaN,NaN];

% Identify baseline period
% CAUTION: 'iter' is used CONTIUOUS for all the following for-loops
for iter = [1:length(transitions)]
    if running_speed(1:transitions(iter))<=moving
        segments.baseline = [1,transitions(iter)];
    else
        break
    end
end

% Identify sample run period
for iter = [iter:length(transitions)]
    % Handle scenarios with and without baseline recordings
    if isnan(segments.baseline)
        segments.sample=[1,transitions(iter)];
    else
        iter=iter-1;
    end

    if any(running_speed(transitions(iter):transitions(iter+1))>moving)
        if isnan(segments.sample)
            segments.sample = [transitions(iter)+1,transitions(iter+1)];
        else
            segments.sample = [segments.sample(1),transitions(iter+1)];
        end
    else
        % Determine direction of sample run based on x-movement direction
        try
            if mean(diff(xcord(segments.sample(1)+2:segments.sample(2)-2)))>0
                segments.sample_dir = 'left'; %CHECK THIS!!!
            elseif mean(diff(xcord(segments.sample(1)+2:segments.sample(2)-2)))<0
                segments.sample_dir = 'right';
            end
        end
        break
    end
end

% Identify delay period between runs
for iter = [iter:length(transitions)]
    if all(running_speed(transitions(iter):transitions(iter+1))<=moving)
        if isnan(segments.delay)
            segments.delay = [transitions(iter)+1,transitions(iter+1)];
        else
            segments.delay = [segments.delay(1),transitions(iter+1)];
        end
    else
        break
    end
end

% Identify match-to-sample run period
for iter = [iter:length(transitions)]
    % Cave: iter+1 may allready go out of bounds
    try
        stoppoint = transitions(iter+1);
    catch
        stoppoint = length(running_speed);
    end
    
    if any(running_speed(transitions(iter):stoppoint)>moving)
        if isnan(segments.match)
            segments.match = [transitions(iter)+1,stoppoint];
        else
            segments.match = [segments.match(1),stoppoint];
        end
    else
        % Determine direction of sample run based on x-movement direction
        if mean(diff(xcord(segments.match(1)+2:segments.match(2)-2)))>0
            segments.match_dir = 'left'; %CHECK THIS!!!
        elseif mean(diff(xcord(segments.match(1)+2:segments.match(2)-2)))<0
            segments.match_dir = 'right';
        end
        break
    end
end

% Post-run, defined as interval beween match-run end to recording end
if ~isnan(segments.match(2))
    if segments.match(2)<length(xcord)
        segments.post = [segments.match(2)+1,length(xcord)];
    end
end

% Classification in match- or non-match trials; Based on trace indicator
try
    xpost = mean(xcord(segments.post(1):segments.post(2)));
    ypost = mean(ycord(segments.post(1):segments.post(2)));
    if (xpost>2.5 & ypost<0.05)
        segments.same_dir = 1;
    elseif (xpost<0.05 & ypost<0.05)
        segments.same_dir = 0;
    end
end
    % Alternative Classification based on sample and match direction.
if ~isfield(segments,'same_dir')    
    try
        if strcmp(segments.sample_dir, segments.match_dir)
            segments.same_dir = 1;
        else
            segments.same_dir = 0;
        end
    end
end
end