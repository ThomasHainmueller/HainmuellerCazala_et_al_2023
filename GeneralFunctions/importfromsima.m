function data = importfromsima(filename)
% Import image data as 5D array from a motion corrected sima file.

%import py.motion_corr.load_simafile
try
    import py.motion_corr.*
catch
    try
        import py.motion_corr.load_simafile.*
    catch
        % Special case for path setting problem
        old = cd('D:/Dropbox/Scanbox_Analysis/');
        import py.motion_corr.load_simafile.*
        cd(old);
    end
end

if exist(strrep(filename,'.sima','_temp.mat'),'file')
    data = load(strrep(filename,'.sima','_temp.mat'),'data');
else
    %py.motion_corr.load_simafile(filename);
    load_simafile(filename);
    data = load(strrep(filename,'.sima','_temp.mat'),'data');
end

data = data.data;

% DEBUGGIN ONLY!
% delete(strrep(filename,'.sima','_temp.mat'))
end