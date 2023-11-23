function data = load_bigsima(it, nframes, savedir)
% Function for importing sima files which are to large for the 'export
% matlab' default option. Accepts a python iterator object (sequence) and
% the number of frames and returns a uint16 array.
% Format: [frame, plane, x, y, channel]

%res = py.motion_corr.create_iterator([olddir,'\',simadir]);
%nframes = cellfun(@int64,cell(res(2)));
%it = cell(res(1)){1}; %Iterator object for sima sequence

for n = 1:nframes
    thisframe = it.next();
    data_size = cellfun(@int64,cell(thisframe.shape));
    thisframe = uint16(py.array.array('d', py.numpy.nditer(...
        thisframe, pyargs('order', 'C'))));
    thisframe = reshape(thisframe,fliplr(data_size));
    thisframe = permute(thisframe,[length(data_size):-1:1]);
    
    if ~exist('data','var')
        data = zeros([nframes,size(thisframe)],'uint16');
    end
    
    data(n,:,:,:,:) = thisframe;
end

if nargin>2
    save(savedir,'data','-v7.3');
end

end