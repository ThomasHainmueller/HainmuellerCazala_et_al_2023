function ephysdata = read_ephysfile(filename,nchannels)
% Extract data from a .ephys file. First row is the number of the present
% image frame for each datapoint. Following rows are the recorded channels.
% 
% nchannels = number of recording channels defined in 'scanbox.mat' PLUS
% ONE (for framenumbers!);

if nargin<2
    nchannels = 'auto';
end

fID = fopen(filename);
ephysdata = fread(fID,'single');
fclose(fID);

% Automatically determine number of channels assuming that the frames (1st)
% channel has a constant value for at least the 10 first frames.
if strcmp(nchannels,'auto') 
    nchan = 1;
    while nchan<length(ephysdata)
        if all(ephysdata(1:nchan:10*nchan) == ephysdata(1))
            nchannels = nchan;
            break
        else
            nchan = nchan+1;
        end
    end
end

ephysdata=reshape(ephysdata,nchannels,length(ephysdata)/nchannels);

end