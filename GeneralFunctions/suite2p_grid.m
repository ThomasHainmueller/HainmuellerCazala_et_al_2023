
saving_name = 'xxx';

%% Run the routine

sbxcategorize();

Bins = 0.1:0.025:2.1; % 80 bins
files = ls('*.ephys');
cell_ID = find(iscell(:,1) == 1); % manually curated cells/boutons
info = load(strrep(files(1,:),'.ephys',''),'info');info = info.info;

if isempty(info.otparam) 
   nplanes = 1;
else
   nplanes = info.otparam(3);
end

if(info.scanmode==0)
   RecordsPerBuffer = info.recordsPerBuffer*2;
else
   RecordsPerBuffer = info.recordsPerBuffer;
end        
 
switch info.channels
       case 1
            info.nchan = 2; % both PMT0 & PMT1
            factor = 1;
       case 2
            info.nchan = 1; % only PMT 0
            factor = 2;
       case 3
            info.nchan = 1; % only PMT 1
            factor = 2;
end

Lines = info.sz(2);
framerate = info.resfreq * (2-info.scanmode) / info.config.lines / nplanes; % calculate acquisition rate
tt=1/framerate;

dataset = [];

% identify accumulative frame numebrs
elapsed_frame = zeros(size(files,1),1);

for run = 1:size(files,1)
    filename = strrep(files(run,:),'.ephys','');   
    d = dir([filename '.sbx']);
    Max_idx =  d.bytes/RecordsPerBuffer/Lines*factor/4 -1; 
    elapsed_frame(run) = round((Max_idx-1)/nplanes); % number of frame per plane
end

% extract behavioral information
for run = 1:size(files,1)
   
    filename = strrep(files(run,:),'.ephys','');      
    
    %% Sort the category first    
    thiscat = load(filename,'category');     
    thiscat = thiscat.category;
    
    % find out how many runs have been stored in this category.
    try
        catrun = length(dataset.metadata.categories{thiscat}.ft)+1;
    catch
        catrun = 1;
    end

    dataset.metadata.categories{thiscat}.filename{catrun} = filename;
    
    %% Time related       
    dataset.metadata.categories{thiscat}.nframes{catrun} = elapsed_frame(run);
    dataset.metadata.categories{thiscat}.acquisition_rate{catrun} = framerate;
    dataset.metadata.categories{thiscat}.frameT{catrun} = 1:1:dataset.metadata.categories{thiscat}.nframes{catrun}; 
    dataset.metadata.categories{thiscat}.ft{catrun} = 1*tt:1*tt:dataset.metadata.categories{thiscat}.nframes{catrun}*tt; 
    
    if run == 1
       dataset.metadata.categories{thiscat}.elapsed_ft{catrun} = 1:1:elapsed_frame(run);
    else
       dataset.metadata.categories{thiscat}.elapsed_ft{catrun} = (sum(elapsed_frame(1:run-1))+1):1:(sum(elapsed_frame(1:run-1))+elapsed_frame(run)); 
    end
    
    %% Import behavour data from .sbx file    
    METADATA = get_behaviourdata(strrep(files(run,:),'.ephys',''));    
    % Location/movement related stuff    
    dataset.metadata.categories{thiscat}.x{catrun} = METADATA(:,1)';
    dataset.metadata.categories{thiscat}.y{catrun} = METADATA(:,2)';
    
    speed_ins = zeros(1,length(data.metadata.categories{1,thiscat}.frameT{catrun}));
    speed_ins(1) = 0;    
    speed_ins(2:end) = diff(data.metadata.categories{1,thiscat}.y{catrun})*framerate;
    index = find(speed_ins > 0);
    [~,edges] = histcounts(speed_ins(index), 100);
    speed_threshold = edges(2);    
    Minspeed = speed_threshold/2.1*4; % speed in the real world --> /2.1*4 (m/s)
    dataset.metadata.categories{thiscat}.moving{catrun}= movingmask2(METADATA(:,2), framerate, Minspeed); 
 
    % Behaviour related stuff: pupil area and licking    
    dataset.metadata.categories{thiscat}.licks{catrun} = METADATA(:,3)';
    dataset.metadata.categories{thiscat}.pupil_area{catrun} = METADATA(:,4)';
    
end                

% extract calcium signals from suite2p files
catrun = [];run = [];

for cell = 1:length(cell_ID)   
    
    dataset.cells{cell}.ID = cell_ID(cell);
    
    for run = 1:size(files,1)
   
    filename = strrep(files(run,:),'.ephys','');    
    
    thiscat = load(filename,'category');     
    thiscat = thiscat.category;
    try
        catrun = length(dataset.cells{cell}.categories{thiscat}.F0)+1;
    catch
        catrun = 1;
    end
    
    if isfield(stat{cell_ID(cell)}, 'iplane') == 1
       dataset.cells{cell}.categories{thiscat}.image_plane{catrun} = stat{cell_ID(cell)}.iplane;
    else
       dataset.cells{cell}.categories{thiscat}.image_plane{catrun} = 0;
    end
            
    F_raw = F(cell_ID(cell),dataset.metadata.categories{thiscat}.elapsed_ft{catrun});
    [dFoT,F0,SIGMA] = normtobaseline2(F_raw, 2,5);
           
    dataset.cells{cell}.categories{thiscat}.dFoT{catrun} = dFoT;
    dataset.cells{cell}.categories{thiscat}.F0{catrun} = F0;
    dataset.cells{cell}.categories{thiscat}.baselineSD{catrun} = SIGMA;
    dataset.cells{cell}.categories{thiscat}.spikes{catrun} = spks(cell_ID(cell),dataset.metadata.categories{thiscat}.frameT{catrun});
                     
    while size(dataset.metadata.categories{thiscat}.y{catrun},2)> size(dataset.cells{cell}.categories{thiscat}.dFoT{catrun},2)
          dataset.metadata.categories{thiscat}.y{catrun}(end)= [];
    end   
            
    while size(dataset.metadata.categories{thiscat}.y{catrun},2)< size(dataset.cells{cell}.categories{thiscat}.dFoT{catrun},2)
          dataset.cells{cell}.categories{thiscat}.dFoT{catrun}(end)= [];
    end   
            
    while size(dataset.metadata.categories{thiscat}.moving{catrun},2)> size(dataset.cells{cell}.categories{thiscat}.dFoT{catrun},2)
          dataset.metadata.categories{thiscat}.moving{catrun}(end)= [];
    end             
            
    dataset.cells{cell}.categories{thiscat}.dFoY{catrun} = ...
    discretize(dataset.cells{cell}.categories{thiscat}.dFoT{catrun}.*dataset.metadata.categories{thiscat}.moving{catrun}, dataset.metadata.categories{thiscat}.y{catrun},Bins);                  

    end
end

data = standard_workflow2(dataset);
save([saving_name,'_swf_for_grid'],'data');
