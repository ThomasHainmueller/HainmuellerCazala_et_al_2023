function affine_alignment(varargin)
% Use Matlab affine transformation to improve alignment of a prealigned
% dataset. Multiply the transformation matrices of the original with the
% improving tform and store the result with the data. Generate new mean
% tiffs. Can also be used in control point selection mode or to just make
% a new tiffstack from the existing tforms. Propagate: Do selected
% alignment procedure once and propagate the obtained modification of the
% tform to all other selected datasets.

%args=struct('alignment_channel','green','selection','all','ControlPointsel',0,'newtiff',0,'propagate',0); 
args=struct('alignment_channel','green','selection','all','mode','CPsel'); 
% selection, e.g. [1 2 4], channels 'green','red','ccim'
% Modes: CPsel, CPselPropagate, affine, newtiff

% Overwrite default parameters if required
for pair = reshape(varargin,2,[])
    if isfield(args,pair{1})
        args.(pair{1})=pair{2};
    else
        error('Input argument does not exist\n')
    end
end

switch args.alignment_channel
    case 'green'
        chan = 1;
    case 'red'
        chan = 2;
    case 'ccim'
        chan = 3;
end

fnames = dir('*.sbx');

if strcmp(args.selection,'all')
    args.selection = 1:length(fnames);
end

% Configure alignment algorithm
[optimizer, metric] = imregconfig('monomodal');
optimizer.MaximumIterations = 1000;
optimizer.MaximumStepLength = .008;

Ref = load(strrep(fnames(1).name,'.sbx','.mat'),'images');
Ref = Ref.images;

for fi = 1:length(args.selection)
    fn = strrep(fnames(args.selection(fi)).name,'.sbx','.mat');

    tform = load(fn,'tform');
    tform = tform.tform;
    
    Mov = load(fn,'images');
    Mov = Mov.images;
    
    for pl = 1:size(Ref,1)
        
        if fi == 1
            for ch = 1:size(Ref,4)
                imwrite(uint16(squeeze(Ref(pl,:,:,ch))),...
                    sprintf('Mean_plane%d_channel%d.tif',pl,ch),...
                    'TIFF','writemode','append')
            end
            continue
        
        else
            switch args.mode
                case 'newtiff'
                    % Nothing to do here
                    
                case 'CPsel'
                    [movpts, fixpts] = cpselect(uint16(squeeze(Mov(pl,:,:,chan))),...
                        uint16(squeeze(Ref(pl,:,:,chan))),'Wait',true);
                    tform{pl} = fitgeotrans(movpts,fixpts,'affine');
                    tform{pl} = maketform('affine',tform{pl}.T);
               
                case 'CPselPropagate'
                    if ~(fi>2)
                    	[movpts, fixpts] = cpselect(uint16(squeeze(Mov(pl,:,:,chan))),...
                            uint16(squeeze(Ref(pl,:,:,chan))),'Wait',true);
                        tform{pl} = fitgeotrans(movpts,fixpts,'affine');
                        tform{pl} = maketform('affine',tform{pl}.T);
                        
                        % Make differential tform capturing modification
                        % for propagation
                        oldtform = load(fn,'tform');
                        oldtform = oldtform.tform;
                        
                        oldA = oldtform{pl}.tdata.Tinv;
                        diffA{pl} = tform{pl}.tdata.T * oldA;
                    else
                        oldA = tform{pl}.tdata.T;
                        tform{pl} = maketform('affine', oldA * diffA{pl});
                    end
                    
                case 'affine'
                    % Find improved transformation matrix
                    thisMov = imtransform(squeeze(Mov(pl,:,:,chan)),tform{pl},...
                        'XData',[1 size(Ref,3)], 'YData', [1 size(Ref,2)]);

                    newtform{pl} = imregtform(uint16(thisMov), uint16(squeeze(...
                        Ref(pl,:,:,chan))), 'affine', optimizer, metric);

                    tform{pl} = maketform('affine',...
                        newtform{pl}.T * tform{pl}.tdata.T); % Composite of transformations
            end
        end
        % transform all image channels, append to tiff stack
        for ch = 1:size(Ref,4)
            thisMov = imtransform(squeeze(Mov(pl,:,:,ch)), tform{pl},...
                'XData',[1 size(Ref,3)], 'YData', [1 size(Ref,2)]);
            
            imwrite(uint16(thisMov),...
                sprintf('Mean_plane%d_channel%d.tif',pl,ch),...
                'TIFF','writemode','append')
        
        end
    end
    save(fn,'tform','-append');
    fprintf('%s aligned\n',fn);    
end