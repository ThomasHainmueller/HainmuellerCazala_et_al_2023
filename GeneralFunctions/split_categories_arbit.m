function newdata = split_categories_arbit(data,category,bounds)
% Make a new datafile in which one category from an original dataset is
% split into arbitrary segment. Input args: category (of the original DS),
% bounds = [n1,n2,n3,...], last run of the 1st, 2nd, 3rd,... segment. 
% --- Run this before 'standard_workflow'. ---

% Mk new dataset w duplicated versions of selected categories.
newdata = regroup_categories(data, 1:length(bounds), ones(1,length(bounds))*category);
bounds = [0,bounds];

for s = 1:length(bounds)-1  
    % splice the metadata section
    fields=fieldnames(data.metadata.categories{category});
    for fi=1:length(fields)
        f = char(fields{fi});
        thisfi = getfield(data.metadata.categories{category}, f);
        
        if strcmp(f,'PV') || strcmp(f,'rewardzone') || strcmp(f,'rewardpoints')
            continue
        end
        % all other subfields of 'categories' can be spliced!
        newdata.metadata.categories{s}=...
            setfield(newdata.metadata.categories{s}, f, thisfi(bounds(s)+1:bounds(s+1)));
    end

    % splice the cells' categorial data
    for n = 1:length(newdata.cells)
        fields = fieldnames(data.cells{n}.categories{category});
        for fi=1:length(fields)
            f = fields{fi};
            thisfi = getfield(data.cells{n}.categories{category}, f);
            newdata.cells{n}.categories{s}=...
                setfield(newdata.cells{n}.categories{s}, f, thisfi(bounds(s)+1:bounds(s+1)));
%             if size(thisfi,2) == nruns
%                 newdata.cells{n}.categories{2*c-1}=...
%                     setfield(newdata.cells{n}.categories{2*c-1}, f ,thisfi(1:2:nruns));
%             else
%                 newdata.cells{n}.categories{2*c-1}=...
%                     setfield(newdata.cells{n}.categories{2*c-1}, f ,thisfi);
            end
        end
    end
    
% splice the general part of each cell's data
%         fields = fieldnames(newdata.cells{n});
%         for fi = 1:length(fields)
%             f = char(fields(fi));
%             % Retarded way of implementation, but all alternatives bug!
%             if strcmp(f,'transientrate') || strcmp(f,'AUCrate') || strcmp(f,'spatialinfo') ||...
%                     strcmp(f,'spatial_P') || strcmp(f,'Placefield_P')
%                 if c==1
%                     rmfield(newdata.cells{n},f);
%                     thisfi = [];
%                 else
%                     thisfi = getfield(newdata.cells{n},f);
%                 end
%                 oldfi = getfield(data.cells{n},f);
%                 thisfi(end+1:end+2) = oldfi(categories(c));
%                 newdata.cells{n} = setfield(newdata.cells{n},f,thisfi);
%             end
%         end

end

