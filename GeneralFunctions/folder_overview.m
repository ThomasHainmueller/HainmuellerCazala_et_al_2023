function folderstack = folder_overview(varargin)
% Convenience function to create a stack of the mean images from all the
% files specified by extension (e.g. '_green_corr.tif', '.sbx') in a folder
% TODO: Implement .sbx opening routines, currently only working for .tif

args=struct('folder',[],'extension','*green_corr.tif');

% Overwrite default parameters if required
for pair = reshape(varargin,2,[])
    args.(pair{1}) = pair{2};
end

if args.folder
    oldfolder = cd(args.folder);
end

GreenFiles = ls(args.extension);
folderstack={};

% 1) Create handles for all complete datasets in folder
for n = size(GreenFiles,1):-1:1
    folderstack{n}.image=mean(ImportMultiTiff(strcat(args.folder,GreenFiles(n,:))),3);
    folderstack{n}.filename=GreenFiles(n,:);
    imwrite(uint16(folderstack{n}.image),'overview.tif','TIFF','writemode','append');
end
end
    
    