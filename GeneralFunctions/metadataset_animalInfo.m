function [animaltab, metads] = metadataset_animalInfo(metads)
% Get the session info (animal ID, number of cells per session) for a
% metadataset in tabularic form. Add on animal ID to metadataset.

for ds = 1:length(metads)
    sessionstr = metads{ds}.metadata.categories{1}.filename{1};
    sessionstr = strsplit(sessionstr,'_');
    restab(ds).date = str2num(sessionstr{1});
    
    if isfield(metads{ds},'animal')
        restab(ds).animal = metads{ds}.animal;
    else
        restab(ds).animal = str2num(sessionstr{2});
        metads{ds}.animal = str2num(sessionstr{2});
    end
    
    restab(ds).ncells = length(metads{ds}.cells);
end

% Convert table to animals
animals = unique([restab.animal]);

for a = 1:length(animals)
    idcs = find([restab.animal] == animals(a));
    animaltab(a).animal = restab(idcs(1)).animal;
    animaltab(a).nsessions = length(idcs);
    animaltab(a).ncells = sum([restab(idcs).ncells]);
end

end