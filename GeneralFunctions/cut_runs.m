function dataset = cut_runs(data)
% Take a dataset where potentially multiple runs along a linear track are
% present and cut the dataset to individual runs. Relies on the
% getlineartracks function.

dataset={};

for n = 1:length(data)
    thisdata = data{n}.data;
    thisfilename = data{n}.filename;
    
    x = thisdata(2,:);
    y = thisdata(3,:);
    
    dx = diff(x);
    dy2 = diff(y).*diff(y);
    
    %find points of aprupt y cor changes (run startpoints)
    rawcutpoints=[find(dy2>0.5),length(y)];
    
    %Iterate trough cutpoints to get periods of at least 50 frames
    cutpoints = [1,rawcutpoints(1)];
    for m = 1:length(rawcutpoints)
        if cutpoints(end)-cutpoints(end-1)<50
            cutpoints(end)=rawcutpoints(m);
        else
            cutpoints(end+1)=rawcutpoints(m);
        end
    end
    
    for o = 1:length(cutpoints)-1
        dataset{end+1}.data = thisdata(:,cutpoints(o):cutpoints(o+1));
        dataset{length(dataset)}.filename = thisfilename;
        
        thisdx = dx(cutpoints(o)+5:cutpoints(o+1)-5);
        %TODO: The run direction currently overrides the category attribute
        if sum(thisdx)>0
            dataset{length(dataset)}.category = 1;
        else
            dataset{length(dataset)}.category = 2;
        end
    end
end
end