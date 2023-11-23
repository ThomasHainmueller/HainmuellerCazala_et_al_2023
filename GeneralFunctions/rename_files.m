function rename_files(fix,varargin)
% Add the string 'fix' to all files in one folder, e.g. for putting the
% date or other information.

args=struct('folder',[],'position','before');

% Overwrite default parameters if required
for pair = reshape(varargin,2,[])
    if isfield(args,pair{1})
        args.(pair{1})=pair{2};
    else
        error('Input argument does not exist\n')
    end
end

files = dir();

for id = 1:length(files)
    if ~files(id).isdir
        f=files(id).name;
        if strcmp(args.position, 'before')
            f = sprintf([fix,f]);
        else
            break
        end
        movefile(files(id).name, f);
    end
end
end