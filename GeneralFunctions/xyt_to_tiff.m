function xyt_to_tiff(f,indices)

if nargin < 2
    % get scans from all subfiles in the mes file.
    indices = [1:length(f)];
end

for index = indices
    [green red metadata] = get_xyt(f,index);
    frametimes = metadata{1};
    xcord = metadata{2};
    ycord = metadata{3};
    green = linewise_bg_sub(green);
    red = linewise_bg_sub(red);
    
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
    xname = strcat(savepath,filename,'_xcord.csv');
    xrawname= strcat(savepath,filename,'_xraw.csv');
    yname = strcat(savepath,filename,'_ycord.csv');
    yrawname = strcat(savepath,filename,'_yraw.csv');
    fname = strcat(savepath,filename,'_frametimes.csv');
    
    % write tif files
    for frame = 1:size(green,3)
       imwrite(green(:,:,frame),gname,'WriteMode','append');
    end
    
    for frame = 1:size(red,3)
        imwrite(red(:,:,frame),rname,'WriteMode','append');
    end
    
    % Write *.csv files for x,y coords and framtimes
    csvwrite(fname,frametimes);
    csvwrite(xname,xcord);
    csvwrite(yname,ycord);
    clear green red frametimes xcord ycord;
end
return
end