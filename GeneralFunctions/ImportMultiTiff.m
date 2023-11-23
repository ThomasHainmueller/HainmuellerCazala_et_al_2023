function[array] = ImportMultiTiff(fname)

info = imfinfo(fname);
num_images = numel(info);

array = imread(fname,1);
for k=2:num_images
    array(:,:,end+1)=imread(fname,k);
end

end
