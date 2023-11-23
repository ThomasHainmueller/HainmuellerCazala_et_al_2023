function extract_folder(folder,subtract_bg)
% Open all .MES files in a given folder and extract their subfiles to TIFF
% plus metadata (intended for Framescans with virtual enviroment data.

if nargin < 2
    subtract_bg = true;
end

oldfolder = cd(folder);
MESfiles = ls('*.mes');


for fileindex = 1:size(MESfiles,1)
    thisfilename = strtrim(MESfiles(fileindex,:));
    success = file_open([],[folder thisfilename]);
    f=mestaghandle('isf');
    
    for subfileindex = 1:length(f)
        try
            s = foldedframe2xyz2(f(subfileindex));
            handleindex = ['f' num2str(subfileindex)];
            s2=BidirectionalZigzagCorr_adaptive_foldedframe(mestaghandle(handleindex), 'pmtUG', {'pmtUG', 'pmtUR'});
            to_tiff(s2(1),1,0,subfileindex,subtract_bg);
            delete(s);
            delete(s2);
            clear s s2;
        catch
            warning([folder thisfilename ' subfile f' num2str(subfileindex) ' could not be exported'])
        end
    end
    clear f;
end
cd(oldfolder);
return
end
