function bg = get_bg_light(folder,varargin)
% Use the deadband (assumed to be at the right-hand image side) to estimate
% the per-frame background light in a recording dataset and save to a
% matlab structure in the target folder.

if nargin<1
    folder = pwd;
end

wdw = 15; % Width of the right-edge deadband in pixels.
verbose = true; % Make a figure of all background fluorescences

sbxfiles = dir('*.sbx');

if verbose
    figure; hold on;
end

for fi = 1:length(sbxfiles)
    try 
        % If the entire file fits into RAM
        vid = sbximport(strrep(sbxfiles(fi).name,'.sbx',''),1,0); % get green channel
        bg.image{fi} = mean(vid,3);
        bg.trace{fi} = squeeze(mean(vid(:,size(vid,2)-wdw:size(vid,2),:),[1,2]));
    
    catch
        % If file is too large
        sbxread(strrep(sbxfiles(fi).name,'.sbx',''),1,1);
        global info
        
        % Get image
        last = min([200 info.max_idx]);
        vid = sbxread(strrep(sbxfiles(fi).name,'.sbx',''),1,last);
        vid = squeeze(vid(1,:,:,:));
        bg.image{fi} = mean(vid,3);
        clear vid
        
        % Get trace
        for fr = info.max_idx:-1:1
            im = sbxread(strrep(sbxfiles(fi).name,'.sbx',''),fr,1);
            trace(fr+1) = squeeze(mean(im(1,:,size(im,3)-wdw:size(im,3)),[2 3]));
        end
        
        bg.trace{fi} = trace;
    end
    bg.filename{fi} = sbxfiles(fi).name;
    
    if verbose
        thisbg = bg.trace{fi};
        plot((thisbg-mean(thisbg))./(2*std(thisbg))-fi);
        drawnow
    end
end

save(fullfile(folder,filesep,'background.mat'),'bg');

end