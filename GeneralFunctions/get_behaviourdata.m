function data = get_behaviourdata(fname, nplanes, ephysrate, display)
% Get behaviourally relevant parameters from the '.ephys' and '_eye.mat'
% file that goes with a '.sbx' recording. 
%
% fname:    filename without extensions
% data:     4-6 x nframes file depending onto whether eyes and licks were
%           monitored with the file.
if nargin<4
    display = 0;
end

if nargin<3
    scanbox_config;
    ephysrate = sbconfig.ephysRate; % Assuming it wasn't changed.
end

% Get nframes of the imaging file
sbxread(fname,1,1); 
global info
nframes = info.max_idx+1; % There is a 0th frame in the imaging file!

if nargin<2
    if info.volscan
        nplanes = info.otparam(3);
    else
        nplanes = 1;
    end
end

ephysdata = read_ephysfile(strcat(fname,'.ephys'),'auto');

% Get a complete vector of frames
vframes = complete_frames(ephysdata(1,:),nframes);
%figure; plot(vframes)

% Get VR x- and y- position
for f = nframes:-1:1
    data(1,f) = mean(ephysdata(2,vframes==f))/2; % Mean x-position value in that frame.
    data(2,f) = mean(ephysdata(3,vframes==f))/2; % Factor 2 conversion for compatibility w. old setup.
end

% Get lick-rate
if size(ephysdata,1) >= 4 % check if lick detection is included
    data(3,:) = count_licks(ephysdata(4,:), ephysrate, vframes, display); % Not typo! data(1,:) = correct frames
else
    data(3,:) = NaN(size(data(2,:)));
end

% Get Pupil area
if exist(strcat(fname,'_eye.mat'),'file')
    try
        eye = load(strcat(fname,'_eye.mat'),'eye');
        eye = eye.eye;
    catch
        eye = sbxeyemotion(strcat(fname,'_eye.mat'));
    end
    [Area{1:length(eye)}] = deal(eye.Area);
    Area = cell2mat(Area);
    
    % Interpolate missing values where pupil hasn't been tracked acuarately
    if length(Area) < round(.6 * nframes) % Special case for bidirectional files before v3.3
        data(4,:) = interp1(find(~isnan(Area))*2, Area(~isnan(Area)), 1:nframes);
    else
        try
            data(4,:) = interp1(find(~isnan(Area)), Area(~isnan(Area)), 1:nframes);
        catch
            data(4,:) = NaN(size(data(3,:)));
        end
    end
else
    data(4,:) = NaN(size(data(3,:)));
end

% Downsample for multi-plane imaging.
if nplanes > 1
    data = downsample2plane(data, nplanes);
else
    data = data'; %Keep dimension conventions
end

end

function frames = complete_frames(frames,nframes)
% Assure that there is a period for each frame of the imaging file in which
% data are sampled to avoid missmatch. Particularly relevant for the
% 'bidirectional' scanning before v3.3 of scanbox.

if max(frames) < round(.6 * nframes) % Special case for bidirectional files before v3.3
    frames = 2*frames; % double frame spacing to correct for missing frames
end

% Assure all frametimes are present
for f = nframes:-1:2
   if ~any(frames==f)
       start = find(frames<f-1,1,'last');
       stop = find(frames>f,1,'first');
       frames(start+round((stop-start)/2):stop-1) = f; 
   end
end

% Special case missing 1, because 0 period is substantial!
if ~any(frames==1)
    interval = length(find(frames==2));
    stop = find(frames==2,1,'first');
    frames(stop-interval:stop-1)=1;
end

frames(~logical((frames>=0).*(frames<=nframes))) = nan; % remove entries that are in no image!

end

function lickcount = count_licks(rawtrace,fs,frames,display)
% Use a raw analog recording trace from arduino capacitive sensing, extract
% the lick-events and return the rate for each frame. 'frames' is a
% steadily growing vector where x is the datapoint and y is the frame in
% which it was recorded.

rawtrace = remove_50Hz(rawtrace,fs); % Filter 50 HZ noise
[b,a] = butter(4,.05,'low'); % use 4th order butterworth to further reduce hf noise
rawtrace = [0,diff(filter(b,a,rawtrace))];
%thresh = 3*std(rawtrace(600:end)); % Use 2 stdev threshold
thresh = .008; % Fixed threshold seems to apply more faithfully for most DS

binlicks = zeros(1,length(frames)); % frames is a vector with as many sampling pts as the ephys file and the value of the current img frame on y

n=600; % Discard 600 datapoints (ringing)
while n < length(rawtrace) 
    if rawtrace(n)>thresh
        binlicks(n)=1;
        n=n+49; % ca. 50 ms interval for 1kHz data
    else
        n=n+1;
    end
end

if display
    figure; plot(binlicks * .1); % 170815 debug
    hold on; plot(rawtrace);
    drawnow
end

for f = max(frames):-1:1
    lickcount(f) = sum(binlicks(frames == f));
end

end

function resampdata = downsample2plane(data, nplanes)
% Cave: Assumes following order 1=x, 2=y, 3=licks, 4=pupil diameter,
% because special treatment (summation) is applied to licks.

maxidx = floor(size(data,2)/nplanes); % Make sure frames are an integer multiple of planes

for r = size(data,1):-1:1
    if r==3
        resampdata(:,r) = sum(reshape(data(r,1:maxidx*nplanes),nplanes,[],1));
    else
        resampdata(:,r) = mean(reshape(data(r,1:maxidx*nplanes),nplanes,[],1));
    end
end
end