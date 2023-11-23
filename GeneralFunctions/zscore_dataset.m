function data = zscore_dataset(data,varargin)
% Take an existing dataset of tdata type and create zscore traces of dFoY
% information.

args=struct('mode','PVcell','overwrite_dFoT',false);
% mode: 'PC', 'PVcell', 'tracewise', 'pause' (= blank periods as baseline).

% Overwrite default parameters if required
for pair = reshape(varargin,2,[])
    if isfield(args,pair{1})
        args.(pair{1})=pair{2};
    else
        error('Input argument does not exist\n')
    end
end

if strcmp(args.mode,'PC')
    return
elseif strcmp(args.mode,'PVcell')
    % Make zscore. sigma calculation based on the 10% lowest continuous 20s
    % segment's SD.
    for n=1:length(data.cells)
        for c = 1:length(data.cells{n}.categories)
            traces = deal(data.cells{n}.categories{c}.dFoT);
            traces = cat(2,traces{:});
            traces = traces(~isnan(traces));
            % Subtract the 8% value for each 75 frames (usually ~15s)
            traces = remove_fluctuations(traces,75,8);
            % Subtract median to assure that baseline is centered to zero
            traces = traces-median(traces);
            % Note the standard deviation of baseline for each category
            [sigma, baseline]=findsigma(traces,100,10);
            % Actual z-scoring.
            for tr = 1:length(data.cells{n}.categories{c}.dFoT)
                if args.overwrite_dFoT
                    data.cells{n}.categories{c}.dFoT{tr}=...
                        (data.cells{n}.categories{c}.dFoT{tr}-baseline)./sigma;
                else
                    data.cells{n}.categories{c}.zscored{tr}=...
                        (data.cells{n}.categories{c}.dFoT{tr}-baseline)./sigma;
                end
            end
            clear traces sigma
        end
    end
elseif strcmp(args.mode,'tracewise')
    for n=1:length(data.cells)
        for c = 1:length(data.cells{n}.categories)
            for tr = 1:length(data.cells{n}.categories{c}.dFoT)
                baseline=nanmean(data.cells{n}.categories{c}.dFoT{tr});
                sigma=nanstd(data.cells{n}.categories{c}.dFoT{tr});
                data.cells{n}.categories{c}.zscored{tr}=...
                    (data.cells{n}.categories{c}.dFoT{tr}-baseline)./sigma;
            end
            clear sigma
        end
    end
elseif strcmp(args.mode,'pause')
    % Traces merged, sigma calculation based on the pause intervals between
    % tracks.
    
    lowest = 0.05; % Lowest y-value that is considered 'on track'
    highest = 2.05; % Highest y-value that is considered 'on track'
    
    for c = 1:length(data.metadata.categories)
        ytrace=deal(data.metadata.categories{c}.y);
        ytrace=cat(2,ytrace{:});
        %ytrace = ytrace(~isnan(ytrace));
        pauses = ytrace<lowest|ytrace>highest; % Pause intervals between tracks
        
        for n=1:length(data.cells)
            traces = deal(data.cells{n}.categories{c}.dFoT);
            traces = cat(2,traces{:});
            %traces = traces(~isnan(traces));
            % Subtract median to assure that baseline is centered to zero
            baseline=nanmedian(traces(pauses));
            sigma=nanstd(traces(pauses));
            % Actual z-scoring.
            for tr = 1:length(data.cells{n}.categories{c}.dFoT)
                data.cells{n}.categories{c}.zscored{tr}=...
                    (data.cells{n}.categories{c}.dFoT{tr}-baseline)./sigma;
            end
            clear traces sigma
        end
    end
end