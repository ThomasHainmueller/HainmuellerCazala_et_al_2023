function categorize_objects(varargin)
% For determine categories for object-place association expts. Object
% presence is indicated as high-level in the x-coordinate. Enter the
% numbers of the runs immediately before the object was changed.

args=struct('object_changed', [0,15,30,45,60], 'nobjects', 2,...
    'locations', [.5, 1.3; 1.9, 2.7; 2.9, 3.6]);
% runs after which object was replaced. Possible object locations (x
% coordinates).

% Overwrite default parameters if required
for pair = reshape(varargin,2,[])
    if isfield(args,pair{1})
        args.(pair{1})=pair{2};
    else
        error('Input argument does not exist\n')
    end
end

files = dir('*.ephys');

for block = 1:length(args.object_changed)-1
    
    for r = args.object_changed(block)+1:args.object_changed(block+1)
        traces = read_ephysfile(files(r).name,'auto');
        
        for loc = 1:size(args.locations,1);
            location = traces(3,:)>args.locations(loc,1) & traces(3,:)<args.locations(loc,2);
            locmean(loc) = mean(traces(2,location));
        end
        
        thispos = find(locmean==max(locmean));
        category = thispos + mod((block-1),args.nobjects)*size(args.locations,1);
        save(strrep(files(r).name,'.ephys','.mat'), 'category', '-append');
    end
end

end

% Old manual function, deprecated 180326
% for n = 1:length(files)
%     fprintf([files(n).name,'\n'])
%     traces = read_ephysfile(files(n).name,'auto');
%     figure; plot(traces(3,:), traces(2,:));
%     category = input('Enter category: \n');
%     save(strrep(files(n).name,'.sbx','.mat'), 'category', '-append');
%     close all
% end