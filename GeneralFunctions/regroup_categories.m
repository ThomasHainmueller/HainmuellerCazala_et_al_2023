function data = regroup_categories(data,newcats,oldcats)
% Resort the categories of an existing tdata type dataset. Insert empty
% 'spacer' categories for 0 values in [newcats] to make datasets compatible
% for merging in metadatasets.

olddata = data; % Make copy of original.
emptycats=newcats(oldcats==0);
newcats=newcats(oldcats~=0); % Leave only non-zero indices for search
oldcats=oldcats(oldcats~=0);

% Regroup traces
for n=1:length(data.cells)
    % Regroup cell parameters
    fields=fieldnames(data.cells{n});
    for fi=1:length(fields)
        f=fields{fi};
        if strcmp(f,'categories')
            data.cells{n}=rmfield(data.cells{n},'categories'); % Make sure only selected data get copied
            data.cells{n}.categories(newcats)=olddata.cells{n}.categories(oldcats);
            for c=emptycats
                data.cells{n}.categories{c}=[];
                data.cells{n}.categories{c}.dFoY{1} = NaN(80,1);
            end
        elseif size(getfield(data.cells{n},f)) == [1,length(oldcats)]
            thisfield=getfield(data.cells{n},f);
            newfield(newcats)=thisfield(oldcats);
            newfield(emptycats)=NaN;
            data.cells{n}=setfield(data.cells{n},f,newfield);
        elseif size(getfield(data.cells{n},f)) == [length(oldcats),length(oldcats)]
            thisfield=getfield(data.cells{n},f);
            newfield(newcats,newcats)=thisfield(oldcats,oldcats);
            newfield(emptycats,:)=NaN;
            newfield(:,emptycats)=NaN;
            data.cells{n}=setfield(data.cells{n},f,newfield);
        end
    end    
end


% Regroup metadata
fields=fieldnames(data.metadata);
for fi=1:length(fields)
    f=fields{fi};
    if strcmp(f,'categories')
        data.metadata=rmfield(data.metadata,'categories');
        data.metadata.categories(newcats)=olddata.metadata.categories(oldcats);
        for c=emptycats
            data.metadata.categories{c}=[];
        end
    elseif size(getfield(data.metadata,f)) == [length(oldcats),length(oldcats)]
        thisfield=getfield(data.metadata,f);
        newfield(newcats,newcats)=thisfield(oldcats,oldcats);
        newfield(emptycats,:)=NaN;
        newfield(:,emptycats)=NaN;
        data.metadata=setfield(data.metadata,f,newfield);    
    elseif size(getfield(data.metadata,f)) == [1,length(oldcats)]
        thisfield=getfield(data.metadata,f);
        newfield(newcats)=thisfield(oldcats);
        newfield(emptycats)=NaN;
        data.metadata=setfield(data.metadata,f,newfield);

    end
end
end
