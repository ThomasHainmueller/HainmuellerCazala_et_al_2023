function [frametimes,x,y] = get_vrposition(filename,nchannels,fsamp,alignmethod,scanmode)
% Extract x- and y- position from a .ephys file. Average position and
% acquisition times (in ms from 1st frame) of each frame.
% Input arguments: nchannels analog input chans +1 for frame number; fsamp
% from the scanbox.config file. Default 1kHz

if nargin < 5
    scanmode = 'unidirectional';
end

if nargin < 4
    alignmethod = 'nonrigid';
end

fID = fopen(filename);
data = fread(fID,'single');
fclose(fID);

data=reshape(data,nchannels,length(data)/nchannels);
if strcmp(alignmethod,'nonrigid')
    imfilename = strrep(filename,'.ephys','_nonrigid');
else
    imfilename = strrep(filename,'.ephys',''); % e.g. for aligned Tiffs
end

sbxread(imfilename,1,1);
global info
nframes = info.max_idx; % Assure compatibility with image data

if strcmp(scanmode,'unidirectional')
    for f = nframes:-1:1
        if all(data(1,:)~=f) % special handle, in case a frame wasn't recorded.
            lastf = data(1,find(data(1,:)<f,1,'last')); % Use next existing frame before f.
            frametimes(f,1) = find(data(1,:)==lastf,1,'first')*1000/fsamp+(f-lastf); % frametime in ms, add 1 ms for each missing frame
            fo = lastf;
        else
            frametimes(f,1)=find(data(1,:)==f,1,'first')*1000/fsamp; % frametime in ms
            fo = f;
        end
        x(f,1) = mean(data(2,data(1,:)==fo))/2; % Mean x-position value in that frame.
        y(f,1) = mean(data(3,data(1,:)==fo))/2; % Factor 2 conversion for compatibility w. old setup.
    end
    
elseif strcmp(scanmode,'bidirectional')
    for f = ceil(nframes/2):-1:1
        if all(data(1,:)~=f) % special handle, in case a frame wasn't recorded.
            lastf = data(1,find(data(1,:)<f,1,'last')); % Use next existing frame before f.
            indexstart = find(data(1,:)==lastf,1,'first');
            indexend = find(data(1,:)==lastf,1,'last');
            indexmid = ceil(indexstart+(indexend-indexstart)/2);            
        else
            indexstart = find(data(1,:)==f,1,'first');
            indexend = find(data(1,:)==f,1,'last');
            indexmid = ceil(indexstart+(indexend-indexstart)/2);
        end
        % Get data for the "forward-" and "backward-" frame separately
        frametimes(2*f-1,1)=indexstart*1000/fsamp; % frametime in ms
        frametimes(2*f,1)=indexmid*1000/fsamp;
        x(2*f-1,1) = mean(data(2,indexstart:indexmid))/2; % Mean x-position value in that frame.
        x(2*f,1) = mean(data(2,indexmid:indexend))/2;
        y(2*f-1,1) = mean(data(3,indexstart:indexmid))/2; % Factor 2 conversion for compatibility w. old setup.
        y(2*f,1) = mean(data(3,indexmid:indexend))/2; 
    end
end

frametimes=frametimes(1:nframes,1); % Make sure data comply with the imagin file!
x=x(1:nframes,1);
y=y(1:nframes,1);
clear data
end