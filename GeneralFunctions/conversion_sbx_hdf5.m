
files = ls('*.sbx');

for run = 1:size(files,1)
    
    file_name = strrep(files(run,:),'.sbx','');
    sbx2h5(file_name)
    
end
