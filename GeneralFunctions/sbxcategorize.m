function categories = sbxcategorize()
% Revise the categories of the datasets in a given folder before signal
% extraction.

win = 100; % Window for X smoothing (category search)
files = dir('*.sbx');
day = 1;
categories = zeros(length(files),2);

for n = 1:length(files)
    thisX = read_ephysfile(strrep(files(n).name,'.sbx','.ephys'),'auto');
    thisX = thisX(2,:);
    thiscat = max(filter(ones(1,win)/win,1,thisX));
    % Categorize dataset according to x-values (maximum 2.6 for 8-bit,
    [files(n).nativecategory, categories(n,2)] = deal(numcategory(thiscat,10,2.7));
    
    if n>1
        lastdate = datetime(files(n-1).datenum,'ConvertFrom','datenum');
        thisdate = datetime(files(n).datenum,'ConvertFrom','datenum');
        if thisdate.Day > lastdate.Day || thisdate.Month > lastdate.Month
            day = day+1;
        end
    end
    
    [files(n).reldate, categories(n,1)] = deal(day);
end

for n=1:length(files)
    files(n).newcategory = (categories(n,1)-1) * length(unique(categories(:,2)))+...
        find(unique(categories(:,2))==categories(n,2));
    
    fprintf('Run: %1$i - %2$s - category: %3$i \n',...
         n, files(n).name, files(n).newcategory);
end

% Manual inspection of the data
decision = input('Categories o.k.? y/n: ','s');

% Write categories as variable into .mat file of the recording
if strcmp(decision,'y')
    for n=1:length(files)
        category = files(n).newcategory;
        save(strrep(files(n).name,'.sbx','.mat'), 'category', '-append');
    end
    
% Manual correction of the categories as needed
elseif strcmp(decision,'n')
    boundaries = input('Select Boundaries of the recording sessions (last index): ');
    
    boundaries=[0,boundaries];
    
    for session = 2:length(boundaries) % Recording session
        rcategories = input([num2str(boundaries(session-1)+1),' to ',...
            num2str(boundaries(session)),' enter old and new categories [old; new]: ']);
        for n=boundaries(session-1)+1:boundaries(session) %Indices of this session
            for c = 1:size(rcategories,2)
                if files(n).newcategory == rcategories(1,c)
                    [files(n).newcategory, category] = deal(rcategories(2,c));
                    save(strrep(files(n).name,'.sbx','.mat'), 'category', '-append');
                end
            end
        end
    end
    
else
    fprintf('Input must be y or n, aborting...');
end

end