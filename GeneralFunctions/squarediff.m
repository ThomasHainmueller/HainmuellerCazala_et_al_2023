function diff2 = squarediff(image)

for h = size(image,1)
    for w = size(image,2)-1
        diffmatrix(h,w) = (image(h,w)-image(h,w+1))^2;
    end
end
diff2 = mean(mean(diffmatrix));
return
end