function [m,T] = prealign(varargin)
% Use image alignment tools of matlab to align the images in an sbx file
% plane by plane. Uses a predefined ROI for alignment.
% Usefull to facilitate sima alignment with smaller
% increments or even completely replace sima alignment.

args=struct('folder',[],'channel',1); % channels: 1=green, 2=red

% Overwrite default parameters if required
for pair = reshape(varargin,2,[])
    if isfield(args,pair{1})
        args.(pair{1})=pair{2};
    else
        error('Input argument does not exist\n')
    end
end

% switch folder if necessary
if isstr(args.folder)
    folder = strrep(args.folder,'\','/');
    if ~strcmp(folder(end),'/')
        folder(end+1)='/';
    end

    cd(folder);
end

files = ls('*.sbx');

% Infer nplanes. Assumes that all files in folder have the same nplanes
sbxread(strrep(files(1,:),'.sbx',''),1,0);
global info;

if info.volscan == 1
    nplanes = info.otparam(3);
else
    nplanes = 1;
end

vid = sbximport(strrep(files(1,:),'.sbx',''),args.channel,0);

% display each plane, ask for ROI
for p = nplanes:-1:1
    figure; imagesc(mean(vid(:,:,p:nplanes:end),3));
    h = imrect(); 
    pos(p,:) = wait(h);
    pos(p,:) = round(pos(p,:));
end

for fi = size(files,1):-1:1
    vid = sbximport(strrep(files(fi,:),'.sbx',''),args.channel,0);
    
    % Align Files plane by plane
    for p = nplanes:-1:1
        % Extract appropriate ROI from frames belonging to that plane
        align = vid(pos(p,2):pos(p,2)+pos(p,4), pos(p,1):pos(p,1)+pos(p,3), p:nplanes:end);
        [tforms, m] = video_motion_correction(align, 80); % 40 Frames for reference
        
        % Create mean image and fill T-vector with shift coordinates
        %alignvid = double.empty(size(vid,1),size(vid,2),length(tforms),0);
        for fr = 1:length(tforms)
            %alignvid(:,:,fr) = imwarp(vid(:,:,p+(fr-1)*nplanes),tforms{fr},...
                %'OutputView',imref2d(size(vid(:,:,1))));
            T(p+(fr-1)*nplanes,1) = round(tforms{fr}.T(3,2)); %xshift
            T(p+(fr-1)*nplanes,2) = round(tforms{fr}.T(3,1)); %yshift
        end
        %figure; imagesc(mean(alignvid,3)); %DEBUG!!!!!!!
        % Create Mean image
        %m(:,:,p) = mean(alignvid,3);
    end
    save(strrep(files(fi,:),'.sbx','.align'), 'm');
    save(strrep(files(fi,:),'.sbx','.align'), 'T', '-append');    
end
    % store m and t in (w x h x p x fi) matrix
end