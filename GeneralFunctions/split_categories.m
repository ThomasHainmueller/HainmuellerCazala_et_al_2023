function newdata = split_categories(data,categories)
% create a new datafile in which the specified categories are split by
% half. Allows for intra-category comparisons, etc.

% Mk new dataset w duplicated versions of selected categories.
newdata = regroup_categories(data, 1:2*length(categories), sort([categories,categories]));
newdata.metadata=rmfield(newdata.metadata,'Placefield_P');

for c = 1:length(categories)
    
    nruns = length(newdata.metadata.categories{2*c}.ft);
    
    for n = 1:length(newdata.cells)
        % splice the general part of each cell's data
        fields = fieldnames(newdata.cells{n});
        for fi = 1:length(fields)
            f = char(fields(fi));
            % Retarded way of implementation, but all alternatives bug!
            if strcmp(f,'transientrate') || strcmp(f,'AUCrate') || strcmp(f,'spatialinfo') ||...
                    strcmp(f,'spatial_P') || strcmp(f,'Placefield_P')
                if c==1
                    rmfield(newdata.cells{n},f);
                    thisfi = [];
                else
                    thisfi = getfield(newdata.cells{n},f);
                end
                oldfi = getfield(data.cells{n},f);
                thisfi(end+1:end+2) = oldfi(categories(c));
                newdata.cells{n} = setfield(newdata.cells{n},f,thisfi);
            end
        end
        
        % splice the 'categories' part of each cell's data
        fields = fieldnames(newdata.cells{n}.categories{2*c});
        for fi=1:length(fields)
            f = fields{fi};
            thisfi = getfield(newdata.cells{n}.categories{2*c}, f);
            if size(thisfi,2) == nruns
                newdata.cells{n}.categories{2*c-1}=...
                    setfield(newdata.cells{n}.categories{2*c-1}, f ,thisfi(1:2:nruns));
                newdata.cells{n}.categories{2*c}=...
                    setfield(newdata.cells{n}.categories{2*c}, f ,thisfi(2:2:nruns));
            else
                newdata.cells{n}.categories{2*c-1}=...
                    setfield(newdata.cells{n}.categories{2*c-1}, f ,thisfi);
                newdata.cells{n}.categories{2*c}=...
                    setfield(newdata.cells{n}.categories{2*c}, f ,thisfi);
            end
        end
    end
    
    % Regroup metadata
    fields=fieldnames(newdata.metadata.categories{2*c});
    for fi=1:length(fields)
        f = char(fields{fi});
        thisfi = getfield(newdata.metadata.categories{2*c}, f);
        
        if strcmp(f,'PV') || strcmp(f,'rewardzone') || strcmp(f,'rewardpoints')
            continue
        end
        % all other subfields of 'categories' can be spliced!
        newdata.metadata.categories{2*c-1}=...
            setfield(newdata.metadata.categories{2*c-1}, f ,thisfi(1:2:nruns));
        newdata.metadata.categories{2*c}=...
            setfield(newdata.metadata.categories{2*c}, f ,thisfi(2:2:nruns));
    end
    newdata.metadata.Placefield_P(:,2*c-1) = data.metadata.Placefield_P(:,categories(c));
    newdata.metadata.Placefield_P(:,2*c) = data.metadata.Placefield_P(:,categories(c));
end

end

% % Mk new dataset w duplicated versions of selected categories.
% newdata = regroup_categories(data, 1:2*length(categories), sort([categories,categories]));
% 
% for c = 1:length(categories)
%     nruns = length(newdata.metadata.categories{2*c}.ft);
%     for n = 1:length(newdata.cells)
%         %nruns = size(newdata.cells{n}.categories{2*c}.dFoT,2);
%         
%         % splice the 'categories' part of each cell's data
%         fields = fieldnames(newdata.cells{n}.categories{2*c});
%         for fi=1:length(fields)
%             f=fields{fi};
%             thisfi = getfield(data.cells{n}.categories{2*c}, f);
%             if size(thisfi,2) == nruns
%                 newdata.cells{n}.categories{2*c-1}=...
%                     setfield(newdata.cells{n}.categories{2*c-1}, f ,thisfi(1:2:nruns));
%                 newdata.cells{n}.categories{2*c}=...
%                     setfield(newdata.cells{n}.categories{2*c}, f ,thisfi(2:2:nruns));
%             else
%                 newdata.cells{n}.categories{2*c-1}=...
%                     setfield(newdata.cells{n}.categories{2*c-1}, f ,thisfi);
%                 newdata.cells{n}.categories{2*c}=...
%                     setfield(newdata.cells{n}.categories{2*c}, f ,thisfi);
%             end
%         end
%     end
%     
%     % Regroup metadata
%     fields=fieldnames(newdata.metadata.categories{2*c});
%     for fi=1:length(fields)
%         f=fields{fi};
%         thisfi = getfield(newdata.metadata.categories{2*c}, f);
%         
%         % all subfields of 'categories' can be spliced!
%         newdata.metadata.categories{2*c-1}=...
%             setfield(newdata.metadata.categories{2*c-1}, f ,thisfi(1:2:nruns));
%         newdata.metadata.categories{2*c}=...
%             setfield(newdata.metadata.categories{2*c}, f ,thisfi(2:2:nruns));
%     end
% end
% 
% end

