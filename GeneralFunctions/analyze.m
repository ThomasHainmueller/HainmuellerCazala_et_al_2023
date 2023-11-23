function [masterstack cells rois] = analyze(folder)
% Function to get all extracted signals from a folder and put them in a
% matrix containing frametimes, xcor, ycord, and signal of each ROI (this
% order) concatenated in one large array.

oldfolder = cd(folder);
SignalFiles = ls('*green_corr.tif');

% Import data from all files in the specified folder.
for n = size(GreenFiles,1):-1:1
    % strip 'green_corr.tif' from filename
    %filename = GreenFiles(n,1:end-15);
    filename = strrep(GreenFiles(n,:),'green_corr.tif','');
    GreenName = GreenFiles(n,:);
    RedName = strcat(filename,'red_corr.tif');
    FramesName = strcat(filename,'frametimes.csv');
    XName = strcat(filename,'xcord.csv');
    YName = strcat(filename,'ycord.csv');
    
    try
        thisgreen = ImportMultiTiff(strcat(folder,GreenName));
        thisred = mean(ImportMultiTiff(strcat(folder,RedName)),3);
        % BG subraction, added 15.03.15
        thisgreen=BGsubt(thisgreen,0.01);
        thisred=BGsubt(thisred,0.01);
        G{n}=thisgreen;
        R{n}=thisred;
        thismeta{1} = csvread(strcat(folder,FramesName));
        thismeta{2} = csvread(strcat(folder,XName));
        thismeta{3} = csvread(strcat(folder,YName));
        % Identify events (match-delay-sample runs) in the x,y data
        thismeta{4} = FindEvents(thismeta);
        meta{n} = thismeta;
        clear thisgreen thisred thismeta;
    catch
        warning([filename ' couldnt be loaded, incomplete dataset or failed eventdetection']);
    end
end

% Register all image stacks
reg{1} = G{1};
masterstack = R{1};

for n = length(G):-1:2
    if size(G{n},3)~=0
        [reg{n}, ThisRefReg] = RegisterStacksRef(G{n},R{n},G{1},R{1});
        masterstack = cat(3,masterstack,ThisRefReg);
    end
end
clear G;

% Crop image area to the part that is common to all runs, create green
% masterstack (uint8) for ROI identification.
[x,y]=imcrop_coords(masterstack);
masterstack=masterstack(x,y,:);
allgreen=[];

for n = length(reg):-1:1
    if size(reg{n},3)~=0
        thisgreen=reg{n};
        thisgreen=thisgreen(x,y,:);
        % Convert all thisgreen stacks to uint8 within the same boundaries
        if exist('minval')
            thisgreen(thisgreen<=minval)=minval;
            thisgreen(thisgreen>=maxval)=maxval;
            thisgreen=uint8(thisgreen/(maxval/255));
        else
            [thisgreen, minval, maxval]=imconvert_uint8(thisgreen,0.001,0.999);
        end
        allgreen=cat(3,allgreen,thisgreen);
        clear thisgreen;
    end
end

%rois=allgreen;
%cells=1;
%return

% Extract ROIs from masterstack
rois = FindRois(allgreen(:,:,:,1));

% Extract event-related traces for all ROIs
for run = length(reg):-1:1
    if size(reg{run},3)~=0
        try
    thisgreen = reg{run}(:,:,:,1);    
    thisft = meta{run}{1};
    thisx = meta{run}{2};
    thisy = meta{run}{3};
    thisev = meta{run}{4};
   
    for roi = size(rois,3):-1:1
        %cells(roi).trace(1) = struct();
        thistrace = ExtractRoi(thisgreen,rois(:,:,roi));
        thistrace = normalize(thistrace,0.05);
        
        % Extract event-related data from current roi and run.
        if ~isnan(thisev.baseline(1))
            %m = numel(cells(roi).trace)+1;
            baselineinterval = [thisev.baseline(1):thisev.baseline(2)];
            cells(roi).trace(5*run-4).data = thistrace(baselineinterval);
            cells(roi).trace(5*run-4).frametimes = thisft(baselineinterval);
            cells(roi).trace(5*run-4).type = 'baseline';
            cells(roi).trace(5*run-4).sample_dir = thisev.sample_dir;
        end
        if ~isnan(thisev.sample(1))
            %m = numel(cells(roi).trace)+1;
            sampleinterval = [thisev.sample(1):thisev.sample(2)];
            cells(roi).trace(5*run-3).data = thistrace(sampleinterval);
            cells(roi).trace(5*run-3).frametimes = thisft(sampleinterval);
            cells(roi).trace(5*run-3).xcord = thisx(sampleinterval);
            cells(roi).trace(5*run-3).ycord = thisy(sampleinterval);
            cells(roi).trace(5*run-3).type = 'sample';
            cells(roi).trace(5*run-3).sample_dir = thisev.sample_dir;
            cells(roi).trace(5*run-3).same_dir = thisev.same_dir;
        end        
        if ~isnan(thisev.delay(1))
            %m = numel(cells(roi).trace)+1;
            delayinterval = [thisev.delay(1):thisev.delay(2)];
            cells(roi).trace(5*run-2).data = thistrace(delayinterval);
            cells(roi).trace(5*run-2).frametimes = thisft(delayinterval); 
            cells(roi).trace(5*run-2).type = 'delay';
            cells(roi).trace(5*run-2).sample_dir = thisev.sample_dir;
            cells(roi).trace(5*run-2).match_dir = thisev.match_dir;
            cells(roi).trace(5*run-2).same_dir = thisev.same_dir;
        end       
        if ~isnan(thisev.match(1))
            %m = numel(cells(roi).trace)+1;
            matchinterval = [thisev.match(1):thisev.match(2)];
            cells(roi).trace(5*run-1).data = thistrace(matchinterval);
            cells(roi).trace(5*run-1).frametimes = thisft(matchinterval); 
            cells(roi).trace(5*run-1).xcord = thisx(matchinterval);
            cells(roi).trace(5*run-1).ycord = thisy(matchinterval);
            cells(roi).trace(5*run-1).type = 'match';
            cells(roi).trace(5*run-1).sample_dir = thisev.sample_dir;
            cells(roi).trace(5*run-1).match_dir = thisev.match_dir;
            cells(roi).trace(5*run-1).same_dir = thisev.same_dir;
        end
        if ~isnan(thisev.post(1))
            %m = numel(cells(roi).trace)+1;
            postinterval = [thisev.post(1):thisev.post(2)];
            cells(roi).trace(5*run).data = thistrace(postinterval);
            cells(roi).trace(5*run).frametimes = thisft(postinterval); 
            cells(roi).trace(5*run).type = 'post';
            cells(roi).trace(5*run).match_dir = thisev.match_dir;
            cells(roi).trace(5*run).same_dir = thisev.same_dir;
        end      
    end
        catch
            warning('One file not added')
        end
    end
end 
end