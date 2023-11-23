function motion_correct_redny(varargin)
%MOTION_CORRECT 
% Correcting motion in dataset using sima (python) as a core routine.
% Aligning the datasets and store the required transformations with the mat
% file that goes with DS. Steps can be repeated separately: 1) Deadspace
% for bidirectional scanning, 2) Sima-motion correction, 3) ccim for ROI
% detection, 4) Dataset alignment, 5) Writing tiff of mean images

args=struct('max_disp',[30,30],'states_retained',50,'correction_channel','red',...
    'repeat_from_step',nan,'maxshift',[100,100],'maxrot',-3:.5:3,...
    'alignment_channel','red','useROI',false,'alignment_mode','XCorr2',...
    'bigfile',false,'mode',1,'nlines',8,'ny',true);
% planes: e.g [1,2,3]; channels: 1=green, 2=red; alignment_mode: XCorr2,
% affine

% Overwrite default parameters if required
for pair = reshape(varargin,2,[])
    if isfield(args,pair{1})
        args.(pair{1})=pair{2};
    else
        error('Input argument does not exist\n')
    end
end

fnames = dir('*.sbx');

for f=1:length(fnames)
    % Step1: Boundaries of deadspace and bidirectional fringes.
    fname=strrep(fnames(f).name,'.sbx','');
    info = load(fname,'info');
    info = info.info;
    if isfield(info,'scanbounds') && ~(args.repeat_from_step<=1)
        continue
    else
        frame = sbxread(fname,1,1);
        if strcmp(args.correction_channel, 'red')
            info.scanbounds=bidirectional_bounds(squeeze(frame(2,:,:)));
        else
            info.scanbounds=bidirectional_bounds(squeeze(frame(1,:,:)));
        end
        save(fname,'info','-append');
    end
end

% Initialize some functions and variables for 'parfor'
meanimages = {};

if ~args.ny
    try
        import py.motion_corr.load_simafile.*
        import py.motion_corr.correct_sbx.*
    end
    try
        import py.motion_corr.*
    end
else
    % This is a specific debug operation for the Computer in NYC
    import py.motion_corr_fixedparam.*
end

parfor f=1:length(fnames)
    fname = fnames(f).name;
    infofi = load(strrep(fname,'.sbx',''));
    
    if isfield(infofi,'images') && ~(args.repeat_from_step<=3) %#ok<PFBNS>
        meanimages{f}=infofi.images;
    else
        if ~isempty(dir(strrep(fname,'.sbx','_temp.mat'))) && ~(args.repeat_from_step<=2)
            % Load temp file if motion correction shouldn't be repeated
            corrseq = load(strrep(fname,'.sbx','_temp.mat'),'data');
            corrseq = corrseq.data;
            fprintf('File %s was loaded from temp file\n',fname)
        else
            % Step2: Create sima files, perform motion correction.
            if ~isempty(dir(strrep(fname,'.sbx','.sima'))) && ~(args.repeat_from_step<=2)
                if ~args.bigfile
                    % Load sima file if motion correction shouldn't be repeated
                    load_simafile(strrep(fname,'.sbx','.sima'));
                    corrseq = load(strrep(fname,'.sbx','_temp.mat'),'data');
                    corrseq = corrseq.data;
                else
                    % Big datasets, load the sequence through py iterator
                    res = load_simafile(strrep(fname,'.sbx','.sima'), args.bigfile);
                    nframes = cellfun(@int64,cell(res(2)));
                    it = cell(res(1)); %Iterator object for sima sequence
                    it = it{1};
                    corrseq = load_bigsima(it, nframes, strrep(fname,'.sbx','_temp.mat'));
                end
                fprintf('File %s was loaded from sima file\n',fname)
            else
                if ~isempty(dir(strrep(fname,'.sbx','.sima'))) && (args.repeat_from_step<=2)
                    rmdir(strrep(fname,'.sbx','.sima'))
                end
                
                if ~args.bigfile
                    %correct_sbx(fname, args.max_disp, args.states_retained,...
                        %args.correction_channel); 170716 mod for 8rows
                    if ~args.ny
                        correct_sbx(fname, args.max_disp, args.states_retained,...
                            args.correction_channel, args.bigfile, args.mode, args.nlines);
                    else
                        % Specific bugfix for the NYC computer
                        rcorrect(fname)
                    end
                    corrseq = load(strrep(fname,'.sbx','_temp.mat'),'data');
                    corrseq = corrseq.data;
                    fprintf('File %s has been motion corrected\n',fname)
                else
                    res = correct_sbx(fname, args.max_disp, args.states_retained,...
                        args.correction_channel, args.bigfile, args.mode, args.nlines); % 170716 modified for 8rows
                    nframes = cellfun(@int64,cell(res(2)));
                    it = cell(res(1)); %Iterator object for sima sequence
                    it = it{1};
                    corrseq = load_bigsima(it, nframes, strrep(fname,'.sbx','_temp.mat'));
                end     
            end
        end
        
        shape = [];
        %pyshape = corrseq.shape;
        
        %for n=1:size(pyshape,2)
            %shape(size(pyshape,2)-n+1)=uint16(py.array.array('L',corrseq.shape(n)))
        %end
        
        %corrseq=uint16(py.array.array('L',py.numpy.nditer(corrseq))); % Convert to uint16 matlab type
        %corrseq=permute(reshape(corrseq,shape),[length(shape):-1:1]); % order: frame, plane, row, column, channel
        
        % Step3: Create ccim, save it as '3rd channel with the mean images in the
        % .mat info file.
        %images = squeeze(mean(corrseq,1)); % order: plane, row, column, channel
        dims = size(corrseq);
        images = reshape(mean(corrseq,1),[dims(2:end),1]);
        nchannels = size(corrseq,5);
        
        for n=1:size(corrseq,2) % iterate over planes
            images(n,:,:,nchannels+1) = uint16(CrossCorrImage(single...
                (permute(squeeze(corrseq(2:end,n,:,:,1)),[2,3,1])))*65535);
        end
        
        meanimages{f}=images;
    end
end

refim = []; % Reference image for alignment.

if strcmp(args.alignment_channel,'green')
    alignchan = 1;
elseif strcmp(args.alignment_channel,'red')
    alignchan = 2;
else
    warning('Alignment channel does not exist!');
    return
end

if strcmp(args.alignment_mode,'affine')
    [optimizer, metric] = imregconfig('monomodal');
    optimizer.MaximumIterations = 1000;
    optimizer.MaximumStepLength = .008;
end

for f = 1:length(fnames)
    fname=strrep(fnames(f).name,'.sbx','');
    images = meanimages{f};
    save(fname,'images','-append');
  
    % Step4: Align the datasets via cross correlation. Easier after parfor
    % because reference file must exist.
    % [~, tforms{m}] = xcorr2_reg(handles{1}.align,handles{m}.align,maxshift,rotation);
    % greenstack(:,:,m)=imtransform(handles{m}.green, tforms{m},...
    %    'XData',[1 size(handles{1}.green,2)],'YData',[1 size(handles{1}.green,1)]);
    infofi = load(fname);
    
    % To disable ROI usage for alignment throughout
    if ~any(args.useROI)
        args.useROI(1:size(images,1)) = 0;
    end
    
    if f == 1
        a = [1,0,0; 0,1,0; 0,0,1];
        % 170424 added planewise tforms
        for pl = size(images,1):-1:1
            tform{pl} = maketform('affine',a);
            
            % Select ROI for alignment
            if args.useROI(pl)
                figure; imagesc(squeeze(images(pl,:,:,alignchan)));
                h = imrect();
                pos(pl,:) = wait(h);
                pos(pl,:) = round(pos(pl,:));
                
                % Use only ROI part from reference
                refim{pl} = squeeze(images(pl,...
                    pos(pl,2):pos(pl,2)+pos(pl,4),...
                    pos(pl,1):pos(pl,1)+pos(pl,3),...
                    alignchan));
            else
                % Use entire first mean image as reference
                refim{pl} = squeeze(images(pl,:,:,alignchan));
            end
            corrimages = images;
        end
    %elseif isfield(infofi,'tform') && ~(args.repeat_from_step<=4)
        %tform = infofi.tform;
        %images = permute(imtransform(permute(images,[2,3,1,4]),tform,...
            %'XData',[1 size(refim,2)], 'YData', [1 size(refim,1)]),[3,1,2,4]);
    % This case should actually never apply as either images AND tforms
    % exist or a new set of both is generated as repeatstep is <= 4
    
    else  
        for pl = size(images,1):-1:1
            if args.useROI(pl)
                thismoving = squeeze(images(pl,...
                    pos(pl,2):pos(pl,2)+pos(pl,4),...
                    pos(pl,1):pos(pl,1)+pos(pl,3),...
                    alignchan));
            else
                thismoving = squeeze(images(pl,:,:,alignchan));
            end
            
            % Inserted 170809
            if strcmp(args.alignment_mode,'XCorr2')
                [~, tform{pl}] = xcorr2_reg(refim{pl}, thismoving,...
                    args.maxshift, args.maxrot);
            elseif strcmp(args.alignment_mode,'affine')
                tform{pl} = imregtform(refim{pl}, thismoving, 'affine',...
                    optimizer, metric);
                tform{pl} = maketform('affine',tform{pl}.T);
            else
                warning('Invalid alignment mode, aborting')
                return
            end
            
            corrimages(pl,:,:,:) = permute(imtransform(permute(images(pl,:,:,:),[2,3,1,4]),tform{pl},...
                'XData',[1 size(corrimages,3)], 'YData', [1 size(corrimages,2)]),[3,1,2,4]);
            %'XData',[1 size(refim{pl},2)], 'YData', [1 size(refim{pl},1)]),[3,1,2,4]);
        end
        fprintf('File %s has been aligned\n',fname)
    end
    
    save(fname,'tform','-append')
    
    % Step5: Write Mean images into tiff files.
    for pl = 1:size(meanimages{f},1)
        for chan = 1:size(meanimages{f},4)
            imwrite(uint16(squeeze(corrimages(pl,:,:,chan))),...
                sprintf('Mean_plane%d_channel%d.tif',pl,chan),'TIFF','writemode','append')
        end
    end
    clear images
end

end

