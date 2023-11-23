function data = update_dataset(data)
% Resolve issues with previous versions of imaging datasets and assure that
% all vectors (x,y,movinyg,transientmask,dFoT are oriented as 1 x n vectors
% and all single-number entries in the dataset are double and not cell.

% Correct metadata
for ca = 1:length(data.metadata.categories)
    for run = 1:length(data.metadata.categories{ca}.ft)
       
        % Correct vector orientation to [1,n]
        for field = {'ft','x','y','moving','licks','pupil_area'}
            field = field{:};
            if isfield(data.metadata.categories{ca},char(field))
                for e = 1:length(data.metadata.categories{ca}.(field))
                    data.metadata.categories{ca}.(field){e}=...
                         reshape(data.metadata.categories{ca}.(field){e},...
                         [1,length(data.metadata.categories{ca}.(field){e})]);
                    %fprintf('%s has been updated\n',field);
                end
            end
        end
        
        % Convert single-value cell entries to double
        for field = {'acquisition_rate'}
            field = field{:};
            if iscell(data.metadata.categories{ca}.(field))
                data.metadata.categories{ca}.(field) =...
                    cat(2,data.metadata.categories{ca}.(field){:});
            end  
        end
    end    
end

% Correct cell data
for n = 1:length(data.cells)
    
    % Convert single-value cell entries to double
    for field = {'transientrate','AUCrate','spatialinfo','spatial_P'}
        field = field{:};
        if isfield(data.cells{n},char(field))
            if iscell(data.cells{n}.(field))
                data.cells{n}.(field) =...
                    cat(2,data.cells{n}.(field){:});
            end   
        end
    end
    
    % Correct vector orientation to [1,n]
    for ca = 1:length(data.cells{n}.categories)
        for run = 1:length(data.cells{n}.categories{ca}.dFoT)
            for field = {'dFoT','transientmask'}
                field = field{:};
                if isfield(data.cells{n}.categories{ca},char(field))
                    for e = 1:length(data.cells{n}.categories{ca}.(field))
                        data.cells{n}.categories{ca}.(field){e}=...
                             reshape(data.cells{n}.categories{ca}.(field){e},...
                             [1,length(data.cells{n}.categories{ca}.(field){e})]);
                        %fprintf('%s has been updated\n',field);
                    end
                end
            end
            
            % Supplement baseline SD for old datasets if necessary
            if ~isfield(data.cells{n}.categories{ca},'baselineSD')
                data.cells{n}.categories{ca}.baselineSD(...
                    1:length(data.cells{n}.categories{ca}.dFoT))=...
                    data.cells{n}.categories{ca}.baselinesigma;
            end
        end 
    end
end
end