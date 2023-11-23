function categories = categorize(data)
% Revise the categories of a dataset generated with get_lineartracks()
% according to user input.

% Display current categories and filenames.
for n=1:length(data)
    %warning([num2str(n),'-',data{n}.filename,'-',num2str(data{n}.category)]);
    fprintf('Run: %1$i - %2$s - category: %3$i \n',...
        n,data{n}.filename,data{n}.category);
    categories(n)=data{n}.category;
end

boundaries = input('Boundaries of the recording sessions (last index): ');

boundaries=[0,boundaries];

for session = 2:length(boundaries) % Recording session
    rcategories = input([num2str(boundaries(session-1)+1),' to ',...
        num2str(boundaries(session)),' enter old and new categories [old; new]: ']);
    for n=boundaries(session-1)+1:boundaries(session) %Indices of this session
        for c = 1:size(rcategories,2)
            if data{n}.category == rcategories(1,c)
                categories(n)=rcategories(2,c);
            end
        end
        data{n}.category=categories(n);
    end
end

% Remove categories with zero values
data(find(categories==0))=[];
        
end