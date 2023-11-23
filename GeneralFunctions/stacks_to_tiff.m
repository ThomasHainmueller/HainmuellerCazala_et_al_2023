function stacks_to_tiff(indices)

f=mestaghandle('isf');

if nargin < 1
    % get scans from all subfiles in the mes file.
    indices = [1:length(f)];
end

for index = indices
    try
        [green red] = get_zstack(f,index);
    
        % Turn images 90 deg so that scanlines are horizontal
        green = uint16(permute(green,[2 1 3]));
        red = uint16(permute(red, [2 1 3]));
      
        % generate and savepath name of the subfile from the .mes file name and location
        savepath = get(f(index),1,'FileName');
        subfolders = strfind(savepath,'\');
    
        % Raw Filename, includes subfolders of the MES file except for the 
        % two topmost (e.g. G:\InVivo\'), '\' are replaced by '_'
        filename = savepath(subfolders(2)+1:end-4);
        filename(filename=='\') = '_';
        filename=strcat(filename,'_f',num2str(index));
        % File path to save file, same as the MES file location
        savepath = savepath(1:subfolders(end));
   
        gname = strcat(savepath,filename,'_green.tif');
        rname = strcat(savepath,filename,'_red.tif');

        % write tif files
        for frame = 1:size(green,3)
            imwrite(green(:,:,frame),gname,'WriteMode','append');
        end
    
        for frame = 1:size(red,3)
            imwrite(red(:,:,frame),rname,'WriteMode','append');
        end

        clear green red
    catch
        warning(['Subfile ',num2str(index),' could not be exported to tiff']);
    end
end
return
end