function mes_to_csv(indices, mesfile)
% Export all AUXin1 channels from a .mes file and store them as .csv in the
% same folder.

if nargin < 2
    mesfile = mestaghandle('isf');
end

if nargin < 1
    indices = 1:length(mesfile);
end

for index = indices
    try
        trace = get(mesfile(index),1,'AUXin1');
        trace = gorconvert(trace);
        
        % generate and savepath name of the subfile from the .mes file name and location
        savepath = get(mesfile(index),1,'FileName');
        subfolders = strfind(savepath,'\');
    
        % Raw Filename, includes subfolders of the MES file except for the 
        % two topmost (e.g. G:\InVivo\'), '\' are replaced by '_'
        filename = savepath(subfolders(2)+1:end-4);
        filename(filename=='\') = '_';
        filename=strcat(filename,'_f',num2str(index));
        % File path to save file, same as the MES file location
        savepath = savepath(1:subfolders(end));
   
        tracename = strcat(savepath,filename,'_ephys.csv');
        csvwrite(tracename,trace);
    catch
        warning(['Subfile ' num2str(index) ' could not be exported. Could possibly contain no ephys data.']);
    end
end

end