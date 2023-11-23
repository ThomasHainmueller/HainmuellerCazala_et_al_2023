function SaveWrite_MultiTiff(array,filename)
% function to make sure that a Multi Page tiff is genereated (at least as
% backup) despite of the shitty windows bug that occasionally denies
% editing of existing files with no warning...

try
    for frame = 1:size(array,3)
       imwrite(array(:,:,frame),filename,'WriteMode','append');
    end
catch
    backupname=strrep(filename,'.tif','_backup.tif');
    try
        for frame = 1:size(array,3)
            imwrite(array(:,:,frame),backupname,'WriteMode','append');
        end        
    catch
        backupnameB=strrep(filename,'.tif','_backup2.tif');
        try
            for frame = 1:size(array,3)
                imwrite(array(:,:,frame),backupnameB,'WriteMode','append');
            end  
        catch
            warning([filename ' could not be saved despite tripple attempt, check folder permissions']);
        end
    end
end
    
return
end