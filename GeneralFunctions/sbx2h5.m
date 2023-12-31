function sbx2h5(fname,varargin)
 
% sbx2h5
% Generates h5 file from sbx files
 
fnh = [fname ,'.h5'];
 
z = sbxread(fname,1,1);
global info;
 
if(nargin>1)
N = min(varargin{1},info.max_idx);
else
N = info.max_idx;
end
 
k = 0;
done = 0;
 
blksize = 200; % block size
 
to_read = min(blksize,N-k);
 
while(~done && to_read>0)
try
q = sbxread(fname,k,to_read);
q = squeeze(q(1,:,:,:)); % extract green channel only
q = permute(q,[2 1 3]);
if(k==0)
h5create(fnh,'/data',[796 512 Inf],'DataType','uint16','ChunkSize',[796 512 to_read]);
h5write(fnh,'/data',q,[1 1 1],[796 512 to_read]);
f = waitbar(0,'Converting to hdf5');
else
h5write(fnh,'/data',q,[1 1 k+1],[796 512 to_read]);
end
catch
done = 1;
delete(f);
end
k = k+to_read;
to_read = min(blksize,N-k);
waitbar(k/N,f);
end
 
delete(f);